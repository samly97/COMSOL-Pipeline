classdef Microstructure
    properties
        id 
        porosity 
        tortuosity
        circles
    end
    methods
        function obj = Microstructure(id, ...
                porosity, tortuosity, circles)
            obj.id = id;
            obj.porosity = num2str(porosity, 3);
            obj.tortuosity = num2str(tortuosity, 3);
            obj.circles = Microstructure.Circles_to_Str(circles);
        end
    end
    methods(Static)
        function circles = Circles_to_Str(circles)
            for c = 1:length(circles)
                circles(c).x = num2str(circles(c).x * 10^6, 3);
                circles(c).y = num2str(circles(c).y * 10^6, 3);
                circles(c).R = num2str(circles(c).R * 10^6, 3);
            end
        end
    end
end