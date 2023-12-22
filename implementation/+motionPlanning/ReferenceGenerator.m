classdef ReferenceGenerator < handle
    %REFERENCEGENERATOR This class generates a nominal reference trajectory for the rest of the pipeline.
    
    properties
        context; % The context for which the reference is generated.

        % Internal
        problem; % contains the settings for the optimization

        % Output
        solutionTumblerFrame; % The solution in the tumbler frame.
        solutionInertialFrame; % The solution in the inertial frame (for debugging).
    end
    
    methods
        function self = ReferenceGenerator(dynamicsContext_)
            %REFERENCEGENERATOR Construct an instance of this class
            %
            % INPUT:
            % context_: The context for which the reference is to be generated.
            %
            % OUTPUT:
            % self: An instance of the ReferenceGenerator class.
            %
            % This constructor initializes the context property with the provided context_,
            % sets up the problem options, defines the problem functions, sets the problem bounds,
            % and computes the final, guess, and solution.
            self.context = dynamicsContext_;

            % set solution method
            self.problem.options(1).nlpOpt = optimoptions('fmincon',...
                'Display','iter',...   % {'iter','final','off'}
                'TolFun',1e-3,...
                'MaxFunEvals',5e4,...
                "Algorithm","interior-point",...
                "ConstraintTolerance",1e-1);   %options for fmincon
            self.problem.options(2).nlpOpt = optimoptions('fmincon',...
                'Display','iter',...   % {'iter','final','off'}
                'TolFun',1e-6,...
                'MaxFunEvals',5e4,...
                "Algorithm","interior-point",...
                "ConstraintTolerance",1e-1);   %options for fmincon
            
            self.problem.options(1).verbose = 2;
            self.problem.options(2).verbose = 2;

            self.problem.options(1).method = 'trapezoid'; % Select the transcription method
            self.problem.options(1).trapezoid.nGrid = 40;  %method-specific options
            withGrad = true;
            self.problem.options(1).nlpOpt.SpecifyConstraintGradient = withGrad;
            self.problem.options(1).nlpOpt.SpecifyObjectiveGradient = withGrad;
            %self.problem.options(1).nlpOpt.CheckGradients = true;
            
            self.problem.options(2).method = 'trapezoid'; % Select the transcription method
            self.problem.options(2).trapezoid.nGrid = 60;  %method-specific options
            self.problem.options(2).nlpOpt.SpecifyConstraintGradient = withGrad;
            self.problem.options(2).nlpOpt.SpecifyObjectiveGradient = withGrad;
            %self.problem.options(2).nlpOpt.CheckGradients = true;


            self.problem.func.dynamics = @(t,x,u) dynamics.odefun(t,x,self.context.parameters,u);
            self.problem.func.pathObj = @(t,x,u)( motionPlanning.pathObjective(u) );
            attitude = precomputeAttitude(self);
            self.context.parameters.attitude = attitude;
            self.problem.func.pathCst = @(t,x,u)( motionPlanning.pathConstraintEllipsoid(t,x,self.context.parameters));
            
            % Problem bounds
            self.problem.bounds.initialTime.low = 0;
            self.problem.bounds.initialTime.upp = 0;
            self.problem.bounds.finalTime.low = self.context.ref.maneuverTime;
            self.problem.bounds.finalTime.upp = self.context.ref.maneuverTime;
            
            state_init = self.context.ref.initialChaserState;
            self.problem.bounds.initialState.low = state_init;
            self.problem.bounds.initialState.upp = state_init;
            

            self.problem.bounds.control.low = -self.context.parameters.u_lim; %-inf;
            self.problem.bounds.control.upp = self.context.parameters.u_lim; %inf;

            computeFinal(self);
            computeGuess(self);
            computeSolution(self);
        end

        function setManeuverDuration(self, maneuverDuration_)
            %SETMANEUVERDURATION Set the maneuver duration
            %
            % INPUT:
            % maneuverDuration_: The new maneuver duration.
            %
            % This method sets the maneuver duration, updates the problem bounds,
            % and recomputes the final, guess, and solution.
            self.context.ref.maneuverTime = maneuverDuration_;
            self.problem.bounds.finalTime.low = maneuverDuration_;
            self.problem.bounds.finalTime.upp = maneuverDuration_;
            computeFinal(self);
            computeGuess(self);
            computeSolution(self);
        end

        function setInitialTumblerState(self, initialTumblerState_)
            %SETINITIALTUMBLERSTATE Set the initial tumbler state
            %
            % INPUT:
            % initialTumblerState_: The new initial tumbler state.
            %
            % This method sets the initial tumbler state and recomputes the final, guess, and solution.
            self.context.ref.initialTumblerState = initialTumblerState_;
            computeFinal(self);
            computeGuess(self);
            computeSolution(self);
        end

        function outputReferenceGenerator = generateDebugReferenceInertialFrame(self, t)
            %GENERATEDEBUGREFERENCEINERTIALFRAME Generate a debug reference in the inertial frame
            %
            % INPUT:
            % t: The time at which to generate the reference.
            %
            % OUTPUT:
            % outputReferenceGenerator: The generated reference.
            %

            N = self.context.mpc.N;
            T = self.context.mpc.T;
            tt = t:T:t+T*N;
            zRefI = interp1(self.solutionInertialFrame.t, self.solutionInertialFrame.z', tt);
            vRefI = interp1(self.solutionInertialFrame.t, self.solutionInertialFrame.v', tt);
            zRefI = zRefI'; vRefI = vRefI';
            valid = ~isnan(vRefI);
            valid = valid(1,:);
            outputReferenceGenerator.time = tt(valid);
            outputReferenceGenerator.xRefI = zRefI(:,valid);
            outputReferenceGenerator.vRefI = vRefI(:,valid);
        end
        
        function outputReferenceGenerator = generateSafeReferenceTumblerFrame(self, t)
            %GENERATESAFEREFERENCETUMBLERFRAME Generate a safe reference in the tumbler frame
            %
            % INPUT:
            % t: The time at which to generate the reference.
            %
            % OUTPUT:
            % outputReferenceGenerator: The generated reference.
            %

            N = self.context.mpc.N;
            T = self.context.mpc.T;
            tt = t:T:t+T*N;

            zRefT = interp1(self.solutionTumblerFrame.t, self.solutionTumblerFrame.z', tt);
            vRefT = interp1(self.solutionTumblerFrame.t, self.solutionTumblerFrame.v', tt);
            zRefT = zRefT'; vRefT = vRefT';
            valid = ~isnan(vRefT);
            valid = valid(1,:);
            outputReferenceGenerator.time = tt(valid);
            outputReferenceGenerator.xRefT = zRefT(:,valid);
            outputReferenceGenerator.vRefT = vRefT(:,valid);
        end

        function outputRefBound = generateForBound(self,t)
            %GENERATEFORBOUND Generate a reference for the next RPI tube
            %
            % INPUT:
            % t: The time at which to generate the reference.
            %
            % OUTPUT:
            % outputRefBound: The generated reference for the next RPI tube
            %
            
            Ttube = self.context.ref.Ttube;
            tt = linspace(t,t+Ttube,10);

            zRefT = interp1(self.solutionTumblerFrame.t, self.solutionTumblerFrame.z', tt);
            vRefT = interp1(self.solutionTumblerFrame.t, self.solutionTumblerFrame.v', tt);
            zRefT = zRefT'; vRefT = vRefT';
            valid = ~isnan(vRefT);
            valid = valid(1,:);
            outputRefBound.time = tt(valid);
            outputRefBound.xRefT = zRefT(:,valid);
            outputRefBound.vRefT = vRefT(:,valid);
        end


    end


    methods (Access = private)

        function computeFinal(self)
            %COMPUTEFINAL Compute the final state
            %
            % This method computes the final state of the system using the initial tumbler state and the maneuver time.
            % The final state is then used to set the bounds for the final state of the problem.
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);
            [~, s_tumbler] = ode45(@(t,s) dynamics.odetumbler(t,s,self.context.parameters), ...
                [0, self.context.ref.maneuverTime], self.context.ref.initialTumblerState, options);
            s_tumbler = s_tumbler';
            wT_fin = s_tumbler(1:3,end);
            qT_fin = s_tumbler(4:7,end);


            pBerth = self.context.parameters.pBerth;
            R_LT = quat.quat2rotm(qT_fin);
            p_final = R_LT*pBerth;
            w_final_L = R_LT*wT_fin - [0; 0; self.context.parameters.OM];
            v_final = cross(w_final_L, p_final);
            state_final = [p_final; v_final];
            self.problem.bounds.finalState.low = state_final;
            self.problem.bounds.finalState.upp = state_final;
        end

        function computeGuess(self)
            %COMPUTEGUESS Compute an initial guess for the trajectory
            %
            % This method computes an initial guess for the trajectory of
            % the system. The trajectory will be used to initiate the
            % direct collocation to solve the optimal control
            nTime = 100;
            t = linspace(0,self.context.ref.maneuverTime,nTime);
            self.problem.guess.time = t;
            p_init = self.problem.bounds.initialState.low(1:3);
            v_init = self.problem.bounds.initialState.low(4:6);
            p_final = self.problem.bounds.finalState.low(1:3);
            v_final = self.problem.bounds.finalState.low(4:6);
            
            kv = 300;
            spline = motionPlanning.MySpline({p_init, p_final, p_final + kv*v_final/norm(v_final)});
            l = linspace(0,1,nTime);
            %vn = linspace(norm(v_init),norm(v_final),nTime);
            v = @(l) (1-l)*norm(v_init) + l*norm(v_final);
            p_guess = zeros(3,nTime);
            v_guess = zeros(3,nTime);
            
            for i = 1:nTime
                p_guess(:,i) = spline.position(l(i));
                tangent = spline.tangent(l(i));
                v_guess(:,i) = v(l(i))*tangent;
            end

            v_guess(:,end) = v_final;
            self.problem.guess.state = [p_guess; v_guess];
            self.problem.guess.control = zeros(3,nTime);
        end

        function computeSolution(self)
            %COMPUTESOLUTION Compute the solution for the problem
            %
            % This method computes the solution for the problem using the optimTraj function.
            % It interpolates the state and control of the solution to generate a reference in the inertial frame.
            % It then rotates this reference to the tumbler frame.
            % The solutions in both frames are stored in the properties solutionInertialFrame and solutionTumblerFrame.
            soln = optimTraj(self.problem);
            zfun = soln(end).interp.state;
            vfun = soln(end).interp.control;
            tt = linspace(0,self.context.ref.maneuverTime, 1000);
            zI = zfun(tt);
            vI = vfun(tt);
            self.solutionInertialFrame.t = tt;
            self.solutionInertialFrame.z = zI;
            self.solutionInertialFrame.v = vI;

            zRefT = zeros(6,length(tt));
            vRefT = zeros(3,length(tt));
            options = odeset('RelTol',1e-6,'AbsTol',1e-6);
            [~,stateRot] = ode45(@(t,x) dynamics.odetumbler(t,x,self.context.parameters), ...
                            tt,self.context.ref.initialTumblerState, options);
            stateRot = stateRot';

            for i = 1:length(tt)
                q = stateRot(4:7,i);
                zRefT(1:3,i) = quat.rotate(zI(1:3,i),quat.conj(q));
                zRefT(4:6,i) = quat.rotate(zI(4:6,i),quat.conj(q));
                vRefT(:,i) = quat.rotate(vI(:,i),quat.conj(q));
            end

            self.solutionTumblerFrame.t = tt;
            self.solutionTumblerFrame.z = zRefT;
            self.solutionTumblerFrame.v = vRefT;
        end

        function attitude = precomputeAttitude(self)
            %PRECOMPUTEATTITUDE Precompute the attitude for the problem
            %
            % This method precomputes the attitude for the problem based on the problem options.
            % It propagates the tumbler to obtain the end attitude and stores it in a cell array.
            % The cell array is returned as the output of the method.
            N = length(self.problem.options);
            nTrials = zeros(1,N);
            for i = 1:N
                nTrials(i) = self.problem.options(i).trapezoid.nGrid;
            end
            attitude = cell(N,1);
            % 1) propagate tumbler obtain end attitude
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);
            for i = 1:N
                tt = linspace(0, self.context.ref.maneuverTime, nTrials(i));
                [~, s_tumbler] = ode45(@(t,s) dynamics.odetumbler(t,s,self.context.parameters), tt, ...
                    self.context.ref.initialTumblerState, options);
                s_tumbler = s_tumbler';
                qT = s_tumbler(4:7,:);
                T = zeros(3,3,nTrials(i));
                for j = 1:nTrials(i)
                    T(:,:,j) = quat.quat2rotm(qT(:,j));
                end
                attitude{i} = T;
            end
        end


    end
end

