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
            obj.porosity = porosity;
            obj.tortuosity = tortuosity;
            obj.circles = circles;
        end
    end
end