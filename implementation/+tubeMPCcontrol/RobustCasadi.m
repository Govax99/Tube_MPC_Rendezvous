classdef RobustCasadi < handle
    %ROBUSTCASADI This class handles robust control problems (tube-MPC) using CasADi.
    
    properties
        context; % The context for which the robust control problem is solved.
        solver; % The solver used for the control problem.
        args; % The arguments for the solver, initial conditions, constraints
    end
    
    methods
        function self = RobustCasadi(dynamicsContext_, outputRobustSet, tightConstraints)
            %ROBUSTCASADI Construct an instance of this class
            %
            % INPUTS:
            % dynamicsContext_: The dynamics context for which the robust control problem is solved.
            % outputRobustSet: The robust sets for the control problem.
            % tightConstraints: The tightened constraints for the control problem.
            %
            % OUTPUT:
            % self: An instance of the RobustCasadi class.
            %
            % This constructor initializes the context property with the provided dynamics context,
            % and calls the initializeSolver method with the provided output robust set and tightened constraints.
            self.context = dynamicsContext_;
            self.initializeSolver(outputRobustSet, tightConstraints);
        end
        
        function initializeSolver(self, outputRobustSet, tightConstraints)
            import casadi.*
            %INITIALIZESOLVER Initialize the solver for the control problem
            %
            % INPUTS:
            % outputRobustSet: The robust sets for the control problem.
            % tightConstraints: The tightened constraints for the control problem.
            %
            % This method initializes the solver for the control problem based on the system dynamics,
            % the robust sets, and the tightened constraints. It sets up the problem and the solver options.
            A = self.context.mpc.sysd.A;
            B = self.context.mpc.sysd.B;
            nStates = length(A);
            nActions = size(B,2);
            upState = tightConstraints.upState; lwState = -upState;
            upAction = tightConstraints.upAction; lwAction = -upAction;

            Xd = SX.sym('x',nStates);
            Ud = SX.sym('u',nActions);
            rhs = A*Xd + B*Ud;
            fdyn = Function('f',{Xd,Ud},{rhs},{'x','u'},{'x+'});
            
            
            N = self.context.mpc.N;
            U = SX.sym('U',nActions,N);
            P = SX.sym('P',nStates + (N+1)*nStates+N*nActions);
            X = SX.sym('X',nStates,(N+1));
            Xreal = P(1:nStates);
            Xref = P(nStates+1:(N+2)*nStates); Xref = reshape(Xref,[nStates, N+1]);
            Uref = P((N+2)*nStates+1:end); Uref = reshape(Uref,[nActions, N]);
            
            obj = 0; % Objective function
            g = [];  % constraints vector
            
            Q = self.context.mpc.Q;
            R = self.context.mpc.R;
            E = self.context.mpc.P; %end matrix 
            
            X0  = X(:,1); % initial state
            
            % initial condition constraints - > initial state must be contained in the mRPI
            A = outputRobustSet{1}.minRPI.A;
            b = outputRobustSet{1}.minRPI.b;
            gt = A*([X0(1); X0(4)] - [Xreal(1); Xreal(4)]) - b;
            g = [g;gt];
            A = outputRobustSet{2}.minRPI.A;
            b = outputRobustSet{2}.minRPI.b;
            gt = A*([X0(2); X0(5)] - [Xreal(2); Xreal(5)]) - b;
            g = [g;gt];
            A = outputRobustSet{3}.minRPI.A;
            b = outputRobustSet{3}.minRPI.b;
            gt = A*([X0(3); X0(6)] - [Xreal(3); Xreal(6)]) - b;
            g = [g;gt];

            % end condition constraints -> steady state (around last ref.
            % point) must be contained in the MRPI

            Xf  = X(:,end); % initial state
            Xs = Xref(:,end); % steady state

            A = outputRobustSet{1}.MaxRPI.A;
            b = outputRobustSet{1}.MaxRPI.b;
            gt = A*([Xf(1); Xf(4)] - [Xs(1); Xs(4)]) - b;
            g = [g;gt];
            A = outputRobustSet{2}.MaxRPI.A;
            b = outputRobustSet{2}.MaxRPI.b;
            gt = A*([Xf(2); Xf(5)] - [Xs(2); Xs(5)]) - b;
            g = [g;gt];
            A = outputRobustSet{3}.MaxRPI.A;
            b = outputRobustSet{3}.MaxRPI.b;
            gt = A*([Xf(3); Xf(6)] - [Xs(3); Xs(6)]) - b;
            g = [g;gt];

            nConstraintsRPI = size(g,1);
            
            for k = 1:N
                obj = obj+(X(:,k) - Xref(:,k))'*Q*(X(:,k) - Xref(:,k)) + ...
                          (U(:,k)-Uref(:,k))'*R*(U(:,k)-Uref(:,k)) ; % calculate obj
                g = [g; X(:,k+1) - fdyn(X(:,k),U(:,k))]; % constraints from dynamics
            end
            obj = obj + (X(:,N+1) - Xref(:,N+1))'*E*(X(:,N+1) - Xref(:,N+1)); % for now Q -> it should be P
            
            % make the decision variable one column  vector
            OPT_variables = [reshape(X,nStates*(N+1),1);reshape(U,nActions*N,1)];
            
            nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);
            
            opts = struct;
            opts.ipopt.max_iter = 2000;
            opts.ipopt.print_level =0;%0,3
            opts.print_time = 0;
            opts.ipopt.acceptable_tol =1e-8;
            opts.ipopt.acceptable_obj_change_tol = 1e-6;
            
            self.solver = nlpsol('solver', 'ipopt', nlp_prob,opts);
            
            self.args = struct;
            
            self.args.lbg(1:nConstraintsRPI,1) = -inf; % mRPI constraint
            self.args.ubg(1:nConstraintsRPI,1) = 0; % mRPI constraint
            
            self.args.lbg(nConstraintsRPI+1:nConstraintsRPI+nStates*N,1) = 0; % dynamics constraints
            self.args.ubg(nConstraintsRPI+1:nConstraintsRPI+nStates*N,1) = 0;
            
            self.args.lbx = [repmat(lwState,[N+1 1]); repmat(lwAction,[N 1])];
            self.args.ubx = [repmat(upState,[N+1 1]); repmat(upAction,[N 1])];
        end

        function [uLaw, xsol, usol] = solve(self, outputRefGenerator, outputEst)
            %SOLVE Solve the control problem
            %
            % INPUTS:
            % outputRefGenerator: The output of the reference generator for the control problem.
            % outputEst: The output  of the state estimation for the control problem.
            %
            % OUTPUTS:
            % uLaw: The control law for the system.
            % xsol: The nominal MPC state solution for the control problem.
            % usol: The nominal MPC control solution for the control problem.
            %
            % This method solves the control problem given the output reference generator and the output state estimation.
            % It returns the control law, the state solution, and the control solution.
            nStates = length(self.context.mpc.sysd.A);
            nActions = size(self.context.mpc.sysd.B,2);
            N = self.context.mpc.N;
            import casadi.*
            self.args.x0 = [reshape(outputRefGenerator.xRefI,[(N+1)*nStates 1]); reshape(outputRefGenerator.vRefI(:,1:N),[N*nActions 1])];
            self.args.p = [outputEst.chaserState; reshape(outputRefGenerator.xRefI,[(N+1)*nStates 1]); reshape(outputRefGenerator.vRefI(:,1:N),[N*nActions 1])];
        
            sol = self.solver('x0', self.args.x0, 'lbx', self.args.lbx, 'ubx', self.args.ubx,...
                'lbg', self.args.lbg, 'ubg', self.args.ubg,'p',self.args.p);
        
            xsol = full(reshape(sol.x(1:nStates*(N+1)),[nStates, N+1]));
            usol = full(reshape(sol.x(nStates*(N+1)+1:end),[nActions, N]));
            if iscolumn(xsol)
                z = xsol;
                v = usol;
            else
                z = xsol(:,1);
                v = usol(:,1);
            end
            uLaw = @(x) v + self.context.mpc.K*(x - z); % ch
        end
    end
end

