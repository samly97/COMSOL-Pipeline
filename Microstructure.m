classdef Microstructure
    properties
        id % label/iteration the microstructure was created at
        porosity 
        tortuosity % calculated continuum tortuosity
        circles % array of RSA generated circles in the microstructure
    end
    methods
        function obj = Microstructure(id, ...
                porosity, tortuosity, circles)
            % The Microstructure class is intented to be encoded into JSON
            % format to retrieve metadata. It also contains utility
            % function(s) to assist with modeling work.
            obj.id = id;
            obj.porosity = porosity;
            obj.tortuosity = tortuosity;
            obj.circles = circles;
        end
        
        function i_1c = Find_i_1C(self, Cs_max, h_cell)
            % Find_i_1C is a utility function that determines the current
            % density (A/m^2), or 1C-rate, which would theoretically
            % discharge the cell in an hour.
            %
            % For a 2D geometry, the circles are taken to be cylinders in
            % the "z" axis. The capacity is determined based on the
            % summation of the individual cylinders' capacity. Then,
            % current is applied at the current collector (h_cell * z).
            %
            % Cs_max: maximum lithium concentration in molecule (mol/m^3)
            % h_cell: the width of the current collector (um)
            %
            % Returns:
            % i_1c: 1C-rate normalized by the current collector (A/m^2)
            
            F = 96485; % Faraday's constant C/mol
            
            % Summation of all the circles' areas.
            area_arr = arrayfun(@(circle) circle.Area(), self.circles);
            area = sum(area_arr);
            
            h_cell = h_cell * 10^-6;
            
            capacity = F * Cs_max * area; % (C)
            i_1c = capacity / h_cell / 3600; % (Ah/m^2)
        end
        
        function obj = Ready_for_JSON(self)
            % Ready_for_JSON converts the numerical data in self.circles to
            % strings. At this point, the Microstructure class is being
            % encoded into JSON.
            self.porosity = num2str(self.porosity);
            self.tortuosity = num2str(self.tortuosity);
            self.circles = Microstructure.Circles_to_Str(self.circles);
            obj = self;
        end
    end
    methods(Static)
        function circles = Circles_to_Str(circles)
            % Circles_to_Str converts the numerical value of the x, y, and
            % R of each circle object to be ready to be formatted into
            % JSON.
            %
            % circles: array of Circle objects which make up the
            % Microstructure.
            for c = 1:length(circles)
                circles(c).x = num2str(circles(c).x * 10^6);
                circles(c).y = num2str(circles(c).y * 10^6);
                circles(c).R = num2str(circles(c).R * 10^6);
            end
        end
    end
end