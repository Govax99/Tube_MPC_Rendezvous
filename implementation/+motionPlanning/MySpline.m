classdef MySpline
%MySpline This class represents a spline curve.
    %   The MySpline class is used to create and manipulate spline curves. 
    %   It provides methods for constructing the spline, updating it, and 
    %   computing the position, tangent, and normal at any point along the spline.
    
    properties
        knotPoints;      % Knot points of the spline
        controlPoints1;  % First set of control points
        controlPoints2;  % Second set of control points
        Nseg;            % Number of segments in the spline
    end
    
    methods
        function self = MySpline(knotPoints_)
            %MySpline Construct an instance of this class
            %   This constructor takes a set of knot points and creates a spline.
            self.knotPoints = knotPoints_;
            self.Nseg = length(knotPoints_)-1;
            self = updateSpline(self);
        end
        

        function p = position(self,l)
            %position Compute the position on the spline at a given length
            %   This method returns the position on the spline at a given length.
            [P0,P1,P2,P3,t] = getSubSegment(self,l);

            p = (1-t)^3*P0 + 3*(1-t)^2*t*P1 + 3*(1-t)*t^2*P2 + t^3*P3;
        end

        function v = tangent(self,l)
            %tangent Compute the tangent to the spline at a given length
            %   This method returns the tangent to the spline at a given length.
            [P0,P1,P2,P3,t] = getSubSegment(self,l);

            v = -3*(1-t)^2*P0 + 3*(1 - 4*t + 3*t^2)*P1 + 3*(2*t - 3*t^2)*P2 + 3*t^2*P3;
            if norm(v) ~= 0
                v = v/norm(v);
            end
        end

        function n = normal(self,l)
            %normal Compute the normal to the spline at a given length
            %   This method returns the normal to the spline at a given length.
            [P0,P1,P2,P3,t] = getSubSegment(self,l);
            d = -3*(1-t)^2*P0 + 3*(1 - 4*t + 3*t^2)*P1 + 3*(2*t - 3*t^2)*P2 + 3*t^2*P3;
            dd = 6*((1-t)*P0 + (-2+3*t)*P1 + (1-3*t)*P2 + t*P3);
            c = cross(dd,d);
            n = cross(d, c);
            if norm(d) ~= 0 && norm(dd) ~= 0
                n = n/(norm(d)*norm(c));
            end
        end
    end

    methods Access(private)
        function self = updateSpline(self)
            %updateSpline Update the spline based on the current knot points
            %   This method updates the control points of the spline based on the current knot points.
            K = self.knotPoints;
        
            % Generate a, b, c, d for Thomas Algorithm
            A = ones(1, self.Nseg);
            B = 4 * ones(1, self.Nseg);
            C = ones(1, self.Nseg);
            D = cell(1,self.Nseg);
        
            % Left most segment
            A(1) = 0;
            B(1) = 2;
            C(1) = 1;
            D{1} = K{1} + 2 * K{2};
        
            for i = 2:self.Nseg-1
                D{i} = 4 * K{i} + 2 * K{i+1};
            end
        
            % Right segment
            A(self.Nseg) = 2;
            B(self.Nseg) = 7;
            C(self.Nseg) = 0;
            D{self.Nseg} = 8 * K{self.Nseg} + K{self.Nseg+1};
        
            self.controlPoints1 = motionPlanning.ThomasAlgorithm(A, B, C, D);
        
            % Compute ControlPoints2
            for i = 1:self.Nseg-1
                self.controlPoints2{i} = 2 * K{i+1} - self.controlPoints1{i+1};
            end
        
            self.controlPoints2{self.Nseg} = 0.5 * (K{self.Nseg+1} + self.controlPoints1{self.Nseg});
        
        end

        function [P0,P1,P2,P3,t] = getSubSegment(self,l)
            %getSubSegment Get the sub-segment of the spline at a given length
            %   This method returns the control points and parameter of the sub-segment of the spline at a given length.
             if (l == self.Nseg)
	            k = fix(self.Nseg);
	            t = 1.0;
            elseif (l < self.Nseg && l >= 0.0)
                k = fix(l)+1;
                t = mod(l,1);
            end

            P0 = self.knotPoints{k};
            P1 = self.controlPoints1{k};
            P2 = self.controlPoints2{k};
            P3 = self.knotPoints{k+1};
        end
    end 


end

