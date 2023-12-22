classdef Satellite
    %SATELLITE Utility visualization for a simple cubic satellite
    
    properties
        position = 1;
        attitudeQ = 1;
        p_dock = [3/2 0 0];
        % points [8 x 3]
        L = 3;
        points = 3 * [0.5 -0.5 0.5; ...
                      0.5 0.5 0.5; ...
                      0.5 0.5 -0.5; ...
                      0.5 -0.5 -0.5; ...
                      -0.5 -0.5 0.5; ...
                      -0.5 0.5 0.5; ...
                      -0.5 0.5 -0.5; ...
                      -0.5 -0.5 -0.5]';
    end

    methods
        function self = Satellite(config)
            %SATELLITE Construct an instance of this class
            if (isfield(config, "p"))
                self.position = config.p;
            end

            if (isfield(config, "q"))
                self.attitudeQ = config.q;
            end

            if (self.position == 1)
                self.position = zeros(3,size(self.attitudeQ,2));
            end
            if (self.attitudeQ == 1)
                self.attitudeQ = repmat([0; 0; 0; 1],1,size(self.position,2));
            end
        end
        
        function plot(self,k)
            %PLOT display satellite at time step "k"
            p = zeros(8,3);
            for i = 1:size(self.points,2)
                p(i,:) = (self.position(:,k) + quat.quat2rotm(self.attitudeQ(:,k)) * self.points(:,i)).';
            end
            
            l1 = [1 1 1 4 4 5 5 2 2 6 8 3]; % l1, l2 first and second point of a side
            l2 = [2 4 5 3 8 8 6 6 3 7 7 7];
            for i = 1:length(l1)
                plot3([p(l1(i),1), p(l2(i),1)], [p(l1(i),2), p(l2(i),2)], [p(l1(i),3), p(l2(i),3)], 'k')
                hold on;
            end
        end
    end
end

