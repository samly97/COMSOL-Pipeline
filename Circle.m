classdef Circle
    properties
        x % (x,) component of center of circle
        y % (,y) component of center of circle
        R % Radius of circle
    end
    methods
        function obj = Circle(x, y, R)
        % Class constructor
        % x: x-coord (um)
        % y: y-coord (um)
        % R: radius (um)
            obj.x = x;
            obj.y = y;
            obj.R = R;
        end
        
        function bool = Overlaps(self, other)
        % Overlaps checks if other circle overlaps with this one
        % self: the current Circle
        % other: the other Circle
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
        
        function area = Area(self)
        % Area calculates the area of the Circle object
            area = pi * self.R^2;
        end
    end
    
    methods(Static)     
        function [mean, std] = particle_stats(circles) 
        % particle_stats determines particle statistics based on an array
        % of Circle objects passed
        %
        % circles: array of Circle objects
        % 
        % Returns
        % mean: average particle radius (um)
        % std: sample standard deviation (um)
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
        
        function eps = porosity(circles, l_e, h_cell)
        % porosity returns the porosity of the electrode based on an array
        % of Circle objects and the length and width (or height) of the
        % electrode. RSA only approaches the specified porosity within a
        % tolerance, so we need to check the "actual" value
        %
        % circles: array of Circle objects bounded within the electrode
        % l_e: length of the electrode (um)
        % h_cell: width (or height) of the electrode (um)
        %
        % Returns
        % eps: porosity of electrode
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