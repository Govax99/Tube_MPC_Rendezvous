classdef RobustSetGenerator < handle
    %ROBUSTSETGENERATOR This class generates all the robust set needed for
    % the Robust MPC control
    
    properties
        context;
        specifics;
        minRPI;
        maxRPI;
    end
    
    methods
        function self = RobustSetGenerator(context_, specifics)
            %ROBUSTSETGENERATOR Construct an instance of this class
            %
            % INPUTS:
            % context_: The context for which the robust set is generated.
            % specifics: Specific parameters for generating the robust set.
            %
            % OUTPUT:
            % self: An instance of the RobustSetGenerator class.
            %
            % This constructor initializes the context and specifics properties with the provided inputs.
            self.context = context_;
            self.specifics = specifics;
        end

        function setDisturbance(self, Polytope, Zonotope)
            %SETDISTURBANCE Set the disturbance for the system
            %
            % INPUTS:
            % Polytope: The polytope for the disturbance.
            % Zonotope: The zonotope for the disturbance.
            %
            % This method sets the polytope and zonotope for the disturbance in the specifics property.
            self.specifics.unc.zonotope.W = Zonotope;
            self.specifics.unc.polytope.W = Polytope;
        end

        function Fs = incomplete_minRPI(self,s)
            %INCOMPLETE_MINRPI Compute an incomplete minimal robust
            %positively invariant (MRPI) set, equivalent to a low accuracy
            %approximation
            %
            % INPUT:
            % s: The power to which the system matrix is raised, the step
            % of the algorithm at which the computation is stopped.
            %
            % OUTPUT:
            % Fs: The incomplete minRPI set.
            %
            % This method computes an incomplete minRPI set for a given
            % number of steps.
            A = self.specifics.sys.Ak;
            W = self.specifics.unc.polytope.W;


            fr = W.A; g = W.b;
            tmp_alfa = zeros(1,length(g));
            for i = 1:length(g)
                tmp_alfa(i) = W.support((A^s).'*fr(i,:).')/g(i);
            end
            Fs = W;
            for i = 1:s-1
                Fs = Fs + affineMap(W,A^i);
            end
        end
        
        function [Finf, s, alfa] = compute_minRPI(self,eps)
            %COMPUTE_MINRPI Compute the minimal robust positively invariant (mRPI) set
            %
            % INPUT:
            % eps: The precision for the computation of the MRPI set.
            %
            % OUTPUTS:
            % Finf: The mRPI set.
            % s: The final step at which the required precision is reached
            % alfa: final parameter of the algorithm
            %
            % This method computes the mRPI set for a given precision.
            A = self.specifics.sys.Ak;
            W = self.specifics.unc.polytope.W;

            n = length(A);
            fr = W.A; g = W.b;
            tmp_alfa = zeros(1,length(g));
            hwPlus = zeros(1,n); hwMin = zeros(1,n);
            e = eye(n);
            s = 0;
            while true
                % Increment s
                s = s + 1;
            
                % Compute alfa0(s)
                for i = 1:length(g)
                    tmp_alfa(i) = W.support((A^s).'*fr(i,:).')/g(i);
                end
                alfa = max(tmp_alfa);
            
                % Compute M(s)
                for j = 1:n
                    hwPlus(j) = hwPlus(j) + W.support((A^(s-1)).'*e(:,j));
                    hwMin(j) = hwMin(j) + W.support(-(A^(s-1)).'*e(:,j));
                end
                M = max([hwMin, hwPlus]);
                if alfa <= eps/(eps + M)
                    break;
                end
            end
            Fs = W;
            for i = 1:s-1
                Fs = Fs + A^i*W;
            end
            Finf = 1/(1-alfa)*Fs;
            self.minRPI = Finf;
        end

        function Finf = compute_zonotope_minRPI(self,n_steps)
            %COMPUTE_ZONOTOPE_MINRPI Compute the minimal robust positively invariant (mRPI) set using a zonotope
            %
            % INPUT:
            % n_steps: The number of steps for the computation.
            %
            % OUTPUT:
            % Finf: The mRPI set.
            %
            % This method computes the mRPI set for a given number of steps using a zonotope.
            % The computed MRPI set is stored in the minRPI property.
            A = self.specifics.sys.Ak;
            W = self.specifics.unc.zonotope.W;

            [Finf] = RPI_1step(A,W,n_steps);
            self.minRPI = Finf;
        end

        function [Oinf, t, Nt] = compute_MaxRPI(self)
            %COMPUTE_MAXRPI Compute the maximal robust positively invariant (MRPI) set
            %
            % OUTPUTS:
            % Oinf: The MRPI set.
            % t: The number of steps for the computation.
            % Nt: The number of semiplanes composing the set
            %
            % This method computes the MRPI set and the time at which it is computed.
            % The computed MRPI set is stored in the maxRPI property.
            Ak = self.specifics.sys.Ak;
            A = self.specifics.sys.A;
            B = self.specifics.sys.B;
            K = self.specifics.sys.K;
            C = self.specifics.sys.C;
            D = self.specifics.sys.D;
            W = self.specifics.unc.polytope.W;
            Y = self.specifics.Ybound.polytope.Y;
            

            % construct right C,D such that y=C*x + D*w (in input we have
            % instead y=C*x + D*u)
            if isprop(self.context,'disturbance') && isfield(self.context.disturbance,'PHI')
                PHI = self.context.disturbance.PHI;
            else
                PHI = zeros(length(A),W.Dim);
            end
            if isprop(self.context,'disturbance') && isfield(self.context.disturbance,'PSI')
                PSI = self.context.disturbance.PSI;
            else
                PSI = zeros(length(A),W.Dim);
            end
            A = A + B*K;
            B = B*K*PHI + PSI;
            C = C + D*K;
            D = D*K*PHI;

            S = Y.A; r = Y.b; M = length(r);
            t = 0; Nt = 0;
            
            % Step 1: compute r0
            for i = 1:M
                s = S(i,:).';
                r(i) = r(i) - W.support(D.'*s);
            end
            if any(r < 0)
                Oinf = 	Polyhedron.emptySet(length(A));
                return;
            end
            
            H = S*C; g = r;
            while (true)
                % Step 2: compute r_t+1
                for i = 1:M
                    s = S(i,:).';
                    r(i) = r(i) - W.support((C*A^t*B).'*s);
                end
                if any(r < 0)
                    t = t + 1;
                    Oinf = Polyhedron.emptySet(length(A));
                    break;
                end
                
                Ot = Polyhedron(H,g);
                % Step 3: compute H_t+1 and g_t+1
                H = [H; S*C*A^(t+1)];
                g = [g; r];
            
                % Step 4: determine if Ot = Ot+1
                Onew = Polyhedron(H,g);
                Nt = size(H,1);
                if (Onew == Ot)
                    Oinf = Ot;
                    break;
                end
            
                % Step 5: increase time count
                t = t+1;
            end
            self.maxRPI = Oinf;
        end

        
        function outputRobustSet = tightenConstraint(self)
            %TIGHTENCONSTRAINT Tighten the constraint of the system
            %
            % OUTPUT:
            % outputRobustSet: The tightened constraints.
            %
            % This method tightens the constraints of the system by computing the MRPI sets and subtracting them from the state and control bounds.
            % The tightened constraints are returned as the output.

            % Only Zonotope for now
            X = self.specifics.Xbound.polytope.W;
            V = self.specifics.Ubound.polytope.W;
            K = self.specifics.sys.K;
            
            [Oinf, ~, ~] = self.compute_MaxRPI();
            n_steps = 33; % obtained by trial and error
            Finf = utils.zonoToPoly(self.compute_zonotope_minRPI(n_steps));

            outputRobustSet.MaxRPI = Oinf;
            outputRobustSet.minRPI = Finf;
            outputRobustSet.Xtight = X - Finf;
            outputRobustSet.Vtight = V - K*Finf;    
        end
    end
end

