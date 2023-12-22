classdef MonteCarloBound < handle
    %MONTECARLOBOUND This class generates a estimate of the disturbance bounds acting on the system
    
    properties
        context;
    end
    
    methods
        function self = MonteCarloBound(context_)
            %MONTECARLOBOUND Construct an instance of this class
            %
            % INPUT:
            % context_: The context for which the reference is to be generated.
            %
            % OUTPUT:
            % self: An instance of the MonteCarloBound class.
            %
            % This constructor initializes the context property with the provided context_,
            self.context = context_;
        end

        function wBound = disturbanceBound(self,outputReferenceGenerator, outputStateEstimation, mode)
            %DISTURBANCEBOUND Estimates the disturbance bound for a system.
            %
            % This function uses Monte Carlo simulation to find a conservative estimate of
            % the correction due to reference update by estimator. It propagates the trajectory
            % given the nominal initial Tumbler state and computes the disturbance bound based on
            % the specified mode.
            %
            % INPUTS:
            %   self: The object instance.
            %   outputReferenceGenerator: The output of the reference generator.
            %   outputStateEstimation: The output of the state estimation.
            %   mode: The mode for computing the disturbance bound. It can be "max", "sigmasigma", or "n*sigma",
            %         where n is a number. If mode is not provided, it defaults to "max".
            %
            % OUTPUTS:
            %   wBound: The estimated disturbance bound.
            
            if nargin < 4
                mode = "max";
            end

            % Name relevant data
            nState = 6;
            parameters = self.context.parameters;
            nTrials = self.context.unc.nTrials;


            time = outputReferenceGenerator.time;
            N = length(time);
            xRefT = outputReferenceGenerator.xRefT;

            %dalfa = self.context.unc.dalfa;
            %dtheta = self.context.unc.dtheta;
            base = self.context.unc.baseOmega;
            covOmega = base.'*base;
            base = self.context.unc.baseQ;
            covQ = base.'*base;

            omega0 = outputStateEstimation.tumblerState(1:3);
            q0 = outputStateEstimation.tumblerState(4:7);


            % Propagate trajectory given nominal initial Tumbler state
            zI = zeros(nState,N);
            options = odeset('RelTol',1e-6,'AbsTol',1e-6);
            [~,stateRot] = ode45(@(t,x) dynamics.odetumbler(t,x,parameters), ...
                            time, [omega0; q0], options);
            stateRot = stateRot';
            for i = 1:N
                zI(1:3,i) = quat.rotate(xRefT(1:3,i), stateRot(4:7,i));
                zI(4:6,i) = quat.rotate(xRefT(4:6,i), stateRot(4:7,i));
            end
            

            % Use Monte Carlo simulation to find a conservative estimate (max) of
            % the correction due to reference update by estimator
            %if ~isfield(self.context.unc, 'dJ')
            %    self.context.unc.dJ = zeros(3);
            %end
            w = zeros(nState, N, nTrials);
            %h = figure; % debug
            nPlot = nTrials/100;
            %colororder(winter(100)) % Set the color order to use the jet colormap
            for i = 1:nTrials
                omegaUnc = mvnrnd(omega0,covOmega,1)';
                qUnc = mvnrnd(q0,covQ,1)'; %qUnc = quat.quatInCone(q0,dalfa,dtheta);
                param_mod = self.context.monteCarloParameters();
                [~,stateRotUnc] = ode45(@(t,x) dynamics.odetumbler(t,x,param_mod), ...
                            time, [omegaUnc; qUnc], options);
                stateRotUnc = stateRotUnc';
                xI = zeros(6,1);
                xIDebug = zeros(size(zI)); %debug
                for j = 1:N
                    xI(1:3) = quat.rotate(xRefT(1:3,j), stateRotUnc(4:7,j));
                    xI(4:6) = quat.rotate(xRefT(4:6,j), stateRotUnc(4:7,j));
                    xIDebug(:,j) = xI; %debug
                    w(:,j,i) = xI - zI(:,j);
                end
                if mod(i,nPlot) == 0
                    %plot3(xIDebug(1,:),xIDebug(2,:),xIDebug(3,:)) %debug
                    %hold on; %debug
                end
            end
            %l = plot3(zI(1,:),zI(2,:),zI(3,:),'r','LineWidth',3); %debug
            %hold on; %debug
            %grid on; xlabel("x [m]"); ylabel("y [m]"); zlabel("z [m]"); legend(l,"Nominal trajectory x_{0:N}") % debug



            if mode == "max"
                wBound = max(abs(w),[],3); %max for all trials
                %wBound = wBound + wDyn;
            elseif mode == "sigmasigma"
                wBound = std(w,0,[2 3]);
            elseif ~isempty(regexp(mode,"\d?sigma", 'once'))
                n = regexp(mode,"(\d?)sigma", 'once', 'tokens');
                if n == ""
                    n = 1;
                else
                    n = str2double(n);
                end
                wBound = n*std(w,0,3);
            end
        end
    end
end

