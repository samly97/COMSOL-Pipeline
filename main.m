PATH = '/Users/SamLy/Desktop/COMSOL API';
FNAME = 'script_generated.mph';

RESULTS_PATH = '/Users/SamLy/Desktop/COMSOL API/Matlab/';

MIN_EPS = 0.4;
MAX_EPS = 0.6;

% Number of microstructures to generate
NUM_GEN = 2;

F = 96485; % Faraday's Constant C/mol

% Max lithium concentration in NMC particle
Cs_max = 48900; % mol/m^3

% Generate NMC particles to these specs;
min_r = 1;
max_r = 10;
clearance = 0;
l_e = 176;
h_cell = 100;

% Cell parameters
l_sep = 52;
C_so = 980;
C_eo = 1000;

% Initial voltage - sort of important!
% Guess (for positive electrode):
% - Discharge: under Equilibrium potential
% - Charge: over Equilibrium potential
Vo = 4.29; 

% Time: Discharge study settings
duration = 3600; 
interval = 300;

% Pre-assign space for Microstructure. To encode into JSON
structures_to_encode = cell(NUM_GEN, 1);

%%%%%%%%%%%
% GENERIC %
%%%%%%%%%%%
eps = zeros(NUM_GEN, 1);
for i = 1:length(eps)
    eps(i) = rand * (MAX_EPS - MIN_EPS) + MIN_EPS;
end

disp(eps)

% We could probably work with the same model file. Can delete the particles
% and regenerate a set of new ones, then run geometry. Everything else is 
% the same.
model = comsol_fns.setup_model(PATH, FNAME);

for i = 1:length(eps)
    fprintf('Iteration: %d\n', i)
    fprintf('Target porosity: %.2f\n', eps(i));
    
    [circles, model] = comsol_fns.generate_particles(model, ...
        min_r, max_r, clearance, eps(i), l_e, h_cell);
    
    % Calculate 1C-rate
    perim = 0;
    area = 0;
    for j = 1:length(circles)
        perim = perim + 2 * pi * circles(j).R;
        area = area + circles(j).Area();
    end
    capacity = F * Cs_max * area/ (h_cell * 10^-6) / 3600 ;% A h/m^2
    i_1c = capacity;
    
    disp(i_1c)
    
    % On first pass, need to set up base model
    if i == 1
        % i_1C will probably not be a constant number
        % make note to change later
        model = comsol_fns.add_model_parameters(model, ...
            h_cell, l_e, l_sep, C_so, C_eo, i_1c);
        model = comsol_fns.create_geometry(model);
        model = comsol_fns.add_model_variables(model);
        model = comsol_fns.add_materials(model);
        model = comsol_fns.create_mesh(model);
        model = comsol_fns.create_voltage_probe(model);
        
        model = comsol_fns.add_electrochem_pdes(model, Vo);
        model = comsol_fns.add_dilute_transport(model);
        
        model = comsol_fns.add_electrochem_study(model, interval, duration);
        model = comsol_fns.add_tortuosity_study(model);
    else
        % On n >= 2th pass, run geometry to 
        model.component('comp1').geom('geom1').run;
    end
    
    % Probably useful to calculate the porosity or void fraction, either for
    % labelling or whatever down the road
    porosity = Circle.porosity(circles, l_e, h_cell);
    
    % Calculate particle statistics here, may be useful to check later if our
    % program is creating "unique enough" particle configurations later.
    [mean, std] = Circle.particle_stats(circles);
    
    fprintf('porosity: %.2f, mean rad (um): %.2f, std (um): %.2f\n', porosity, ...
        mean*10^6, std*10^6)
    
    model = comsol_fns.run_electrochem_study(model);
    model = comsol_fns.run_tortuosity_study(model);
    
    % Get derived value
    [flux, model] = comsol_fns.flux_for_tortuosity(model);
    tortuosity = porosity * C_eo/(l_e * 10^-6 * flux);
    
    fprintf('tortuosity %.2f\n', tortuosity)

    % Export geometry png here
    model = comsol_fns.export_geometry_pic(model, ...
        sprintf('%s/results/microstructure%d', RESULTS_PATH, i));
    
    % Export electrochemical data here
    model = comsol_fns.export_electrochem_data(model, i, RESULTS_PATH);
    
    % Try removing the particle sequence instead
    model.geom('part1').feature.clear;
    
    % Add Microstructure to data array for encoding
    structures_to_encode{i} = Microstructure(i, porosity, ...
        tortuosity, circles);
end

% Encode to JSON here
json = jsonencode(structures_to_encode);

% Save JSON data to results
fid = fopen(sprintf('%s/results/metadata.json', RESULTS_PATH),'w');
fprintf(fid,'%s\n', json);
fclose(fid);

disp('break here');