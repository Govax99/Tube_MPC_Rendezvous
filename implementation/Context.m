classdef Context < handle
    %CONTEXT This class handles the context for a control problem.
    
    properties
        parameters; % The nominal parameters for the control problem.
        ref; % Reference for the control problem.
        mpc; % Model predictive control for the control problem.
        unc; % Estimated uncertainty for the control problem.
        rpi; % Robust positively invariant sets for the control problem.
        disturbance; % specify matrix PHI PSI for definition of disturbed linear system
        real; % Real parameters for the control problem.
    end
    
    methods
        function param_real = realParameters(self)
            %REALPARAMETERS Get the real parameters for the control problem
            %
            % OUTPUT:
            % param_real: The real parameters for the control problem.
            %
            param_real = self.parameters;
            param_real.J_T = self.real.J_T;
        end

        function param_est = estimatedParameters(self)
            %ESTIMATEDPARAMETERS Get the estimated parameters for the control problem
            %
            % OUTPUT:
            % param_est: The estimated parameters for the control problem.
            %
            param_est = self.parameters;
            param_est.J_T = self.real.J_T + self.unc.dJ.*utils.randnsym(3);
        end

        function param_mod = monteCarloParameters(self)
            %MONTECARLOPARAMETERS Get the parameters for the control
            %problem to use in the Monte Carlo simulation for uncertainty estimation
            %
            % OUTPUT:
            % param_mod: A particular realization of the parameters
            %
            param_mod = self.parameters;
            param_mod.J_T = self.parameters.J_T + self.unc.dJ.*utils.randnsym(3);
        end
    end
end

