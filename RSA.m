function circle_arr = RSA(min_r, max_r, clearance, eps, l_e, h_cell)
% References:
% - https://porespy.org/modules/generated/porespy.generators.RSA.html
% - American Journal of Physics 86, 772 (2018); doi: 10.1119/1.5049954
%
% min_r: minimum electrode particle radius (um)
% max_r: maximum electrode particle radius (um)
% clearance: added distance between particles in case of overlap
% eps: targeted porosity (1 - area_frac_particles/area_cell)
% l_e: length of electrode (um)
% h_cell: height of cell (um)
%
% Returns
% circle_arr: an array of contained Circle objects generated uniformally
% using Random Sequential Addition.
    
    % Convert quantities to (um) scale
    [min_r, max_r, clearance, l_e, h_cell] = to_micro(...
        min_r, max_r, clearance, l_e, h_cell);

    circle_arr = [];
    
    % Keep track of current volume fraction
    vol_frac = 0;
    
    % Area of cell to get volume fraction
    cell_area = l_e * h_cell;
    
    MAX_ITERS = 10000;
    i = 1;
    
    % Arbitrary tolerance, so when porosity is within "x = 0.02", the we
    % stop attempting this algorithm
    while abs(1 - (1 - vol_frac)/eps) > 0.02 && i < MAX_ITERS
        i = i + 1;
        
        R = rand * (max_r - min_r) + min_r;
        x = rand * l_e;
        y = rand * h_cell;
        
        % Check overlap before inserting!
        new_circle = Circle(x, y, R);
        
        % First circle; nothing to compare to
        if isempty(circle_arr)
            circle_arr = [circle_arr new_circle];
            vol_frac = vol_frac + new_circle.Area()/cell_area;
            continue;
        end
        
        % If at any point overlaps, need to get new values, otherwise
        % insert into array
        temp_circle = Circle(x, y, R + clearance);
        overlaps = false;
        for i = 1:length(circle_arr)
            to_compare = circle_arr(i);
            if temp_circle.Overlaps(to_compare)
                overlaps = true;
                break;
            end
        end
        
        % Check if circles contained within the box, so particles don't get
        % cut off when creating geometry.
        if x - R < 0 || x + R > l_e
            continue
        end
        if y - R < 0 || y + R > h_cell
            continue
        end
        
        % Inserting circle here...
        if overlaps == false
            circle_arr = [circle_arr new_circle];
            vol_frac = vol_frac + new_circle.Area()/cell_area;
        end
    end
end

function [min_r, max_r, clearance, l_e, h_cell] = to_micro(...
    min_r, max_r, clearance, l_e, h_cell) 
    min_r = min_r * 1e-6;
    max_r = max_r * 1e-6;
    clearance = clearance * 1e-6;
    l_e = l_e * 1e-6;
    h_cell = h_cell * 1e-6;
end