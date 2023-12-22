classdef Simulator < handle
    %SIMULATORFINAL This class contains the simulator and mock estimation process.
    
    properties
        chaserState; % The state of the chaser in the simulation.
        tumblerState; % The state of the tumbler in the simulation.
        context; % The context of the simulation.
        oldEst; % The previous (most-recent) estimation of the state.
    end
    
    methods
        function self = Simulator(context_)
            %SIMULATORFINAL Construct an instance of this class
            %   This constructor initializes the simulator with a given context.
            self.context = context_;
            self.chaserState = self.context.real.initialChaserState;
            self.tumblerState = self.context.real.initialTumblerState;
            self.reset();
        end

        function outputStateEstimation = reset(self)
            %RESET Resets the state of the simulator.
            %   This function resets the chaser and tumbler states to their initial values.
            
            self.chaserState = self.context.real.initialChaserState;
            self.tumblerState = self.context.real.initialTumblerState;

            % estimation
            [chaserStateEst, tumblerStateEst] = self.estimationProcess();
            outputStateEstimation.parameters = self.context.estimatedParameters();
            outputStateEstimation.chaserState = chaserStateEst; % should be slightly different
            outputStateEstimation.tumblerState = tumblerStateEst;
            self.oldEst = outputStateEstimation;
        end

        function outputStateEstimation = stepUncertain(self,dtEstimation, uLaw)
            %STEPUNCERTAIN Performs a step in the simulation with uncertainty.
            %   This function performs a step in the simulation considering a given control law and time step.

            uAction = uLaw(self.oldEst.chaserState);
            options = odeset('RelTol',1e-12,'AbsTol',1e-12);
            param_real = self.context.realParameters();
            [~, xOut] = ode45(@(t,s) dynamics.odefun(t,s,self.context.parameters,uAction), [0, dtEstimation], self.chaserState, options);
            [~, xOutT] = ode45(@(t,s) dynamics.odetumbler(t,s,param_real), [0, dtEstimation], self.tumblerState, options);

            self.chaserState = xOut(end,:)';
            self.tumblerState = xOutT(end,:)';

            % estimation
            [chaserStateEst, tumblerStateEst] = self.estimationProcess();
            outputStateEstimation.parameters = self.context.estimatedParameters();
            outputStateEstimation.chaserState = chaserStateEst; % should be slightly different
            outputStateEstimation.tumblerState = tumblerStateEst;
            outputStateEstimation.uAction = uAction; % for debug
            self.oldEst = outputStateEstimation;
        end

        function [chaserStateEst, tumblerStateEst] = estimationProcess(self)
            %ESTIMATIONPROCESS Estimates the state of the system (mock version).
            %   This function estimates the state of the chaser and the tumbler based on their current states.
            omega = self.tumblerState(1:3);
            q = self.tumblerState(4:7);
            %dalfa = self.context.unc.dalfa;
            %dtheta = self.context.unc.dtheta;
            base = self.context.unc.baseOmega;
            covOmega = base.'*base;
            base = self.context.unc.baseQ;
            covQ = base.'*base;

            chaserStateEst = self.chaserState + diag([0.001 0.001 0.001 0.00001 0.00001 0.00001])*randn(6,1);
            %qEst = quat.quatInCone(q,dalfa,dtheta);
            qEst = mvnrnd(q,covQ,1)';
            omegaEst = mvnrnd(omega,covOmega,1)';
            tumblerStateEst = [omegaEst; qEst];
        end
    end
end

