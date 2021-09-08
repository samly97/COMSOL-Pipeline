classdef Circle
    properties
        x % (x,) component of center of circle
        y % (,y) component of center of circle
        R % Radius of circle
    end
    methods
        % Class constructor
        function obj = Circle(x, y, R)
            obj.x = x;
            obj.y = y;
            obj.R = R;
        end
        
        % Checks if other circle overlaps with this one
        function bool = Overlaps(self, other)
            dx = self.x - other.x;
            dy = self.y - other.y;
            dist = sqrt(dx^2 + dy^2);
            if dist > self.R + other.R
               bool = false;
               return
            end
            bool = true;
            return
        end
        
        % Area of circle
        function area = Area(self)
            area = pi * self.R^2;
        end
    end
    
    methods(Static)     
        % particle statistics
        % - mean
        % - (sample )std deviation
        function [mean, std] = particle_stats(circles) 
            N = length(circles); 
            
            running_sum = 0;
            for i = 1:N
                running_sum = running_sum + circles(i).R;
            end
            mean = 1/N * running_sum;
            
            sum_squared_diff = 0;
            for i = 1:N
                sum_squared_diff = sum_squared_diff + ...
                    (sum_squared_diff - mean)^2;
            end
            
            std = sqrt(1/(N-1) * sum_squared_diff);
        end
        
        % Porosity of electrode
        function eps = porosity(circles, l_e, h_cell)
            l_e = l_e * 10^-6;
            h_cell = h_cell * 10^-6;
            cell_area = l_e * h_cell;
            area = 0;
            for i = 1:length(circles)
                area = area + pi * circles(i).R^2;
            end
            eps = (cell_area - area)/cell_area;
        end
    end
end