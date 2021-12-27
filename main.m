PATH = '/Users/SamLy/Desktop/COMSOL API';
FNAME = 'script_generated.mph';

RESULTS_PATH = '/Users/SamLy/Desktop/COMSOL API/Matlab';

MIN_EPS = 0.4;
MAX_EPS = 0.6;

% Number of microstructures to generate
NUM_GEN = 5;

% Max lithium concentration in NMC particle
Cs_max = 48900; % mol/m^3

% Generate NMC particles to these specs;
min_r = 1;
max_r = 10;
clearance = 0.2;
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

% Studied C-rates
C_rates = [1/3, 1, 2, 4, 6, 9];

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
    
    % Probably useful to calculate the porosity or void fraction, either for
    % labelling or whatever down the road
    porosity = Circle.porosity(circles, l_e, h_cell);
    
    % Calculate particle statistics here, may be useful to check later if our
    % program is creating "unique enough" particle configurations later.
    [mean, std] = Circle.particle_stats(circles);
    
    fprintf('porosity: %.2f, mean rad (um): %.2f, std (um): %.2f\n', porosity, ...
        mean*10^6, std*10^6)
    
    % Calculate 1C-rate
    microstructure = Microstructure(i, porosity, 0, circles); % update 
    % tortuosity after
    i_1c = microstructure.Find_i_1C(Cs_max, h_cell);
    fprintf('From microstructure: %f\n', i_1c);
    
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
        
        model = comsol_fns.add_electrochem_study(model);
        model = comsol_fns.add_tortuosity_study(model);
    else
        % On n >= 2th pass, run geometry to 
        model.component('comp1').geom('geom1').run;
    end
    
    % For a microstructure, discharge the cell at different C-rates
    for k = 1:length(C_rates)
        C = C_rates(k);
        
        duration = 3600 / C;
        interval = duration / 12;
        
        model = comsol_fns.modify_applied_current(model, i_1c, C);
        model = comsol_fns.setup_electrochem_study_duration(model, ...
            interval, duration);
        model = comsol_fns.run_electrochem_study(model);
        
        % Export electrochemical data here
        model = comsol_fns.export_electrochem_data(model, i, C, RESULTS_PATH);
    end
    
    model = comsol_fns.run_tortuosity_study(model);
    
    % Get derived value
    [flux, model] = comsol_fns.flux_for_tortuosity(model);
    tortuosity = porosity * C_eo/(l_e * 10^-6 * flux);
    
     fprintf('tortuosity %.2f\n', tortuosity)

    % Export geometry png here
    model = comsol_fns.export_geometry_pic(model, ...
        sprintf('%s/results/microstructure%d', RESULTS_PATH, i));
    
    % Try removing the particle sequence instead
    model.geom('part1').feature.clear;
    
    % Update microstructure's tortuosity before encoding to JSON
    microstructure.tortuosity = tortuosity;
    
    % Add Microstructure to data array for encoding
    structures_to_encode{i} = microstructure.Ready_for_JSON();
end

% Encode to JSON here
json = jsonencode(structures_to_encode);

% Save JSON data to results
fid = fopen(sprintf('%s/results/metadata.json', RESULTS_PATH),'w');
fprintf(fid,'%s\n', json);
fclose(fid);

disp('break here');