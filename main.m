PATH = '/Users/SamLy/Desktop/COMSOL API';
FNAME = 'script_generated.mph';

MIN_EPS = 0.4;
MAX_EPS = 0.6;

% Generate NMC particles to these specs;
min_r = 1;
max_r = 10;
clearance = 0;
l_e = 176;
h_cell = 100;

% REMOVE THIS LATER
eps = 0.58;

% Cell parameters
l_sep = 52;
C_so = 980;
C_eo = 1000;
i_1c = 10; % not real number - replace later

% Initial voltage - sort of important!
% Guess (for positive electrode):
% - Discharge: under Equilibrium potential
% - Charge: over Equilibrium potential
Vo = 4.29; 

% Time: Discharge study settings
duration = 30; 
interval = 2;

%%%%%%%%%%%
% GENERIC %
%%%%%%%%%%%
% eps = zeros(5, 1);
% for i = 1:length(eps)
%     eps(i) = rand * (MAX_EPS - MIN_EPS) + MIN_EPS;
% end

% We could probably work with the same model file. Can delete the particles
% and regenerate a set of new ones, then run geometry. Everything else
% should be the same.
%
% But CHECK IT; don't assume
model = comsol_fns.setup_model(PATH, FNAME);

[circles, model] = comsol_fns.generate_particles(model, ...
    min_r, max_r, clearance, eps, l_e, h_cell);

% Probably useful to calculate the porosity or void fraction, either for
% labelling or whatever down the road
porosity = Circle.porosity(circles, l_e, h_cell);

% Calculate particle statistics here, may be useful to check later if our
% program is creating "unique enough" particle configurations later.
[mean, std] = Circle.particle_stats(circles);

model = comsol_fns.add_model_parameters(model, ...
     h_cell, l_e, l_sep, C_so, C_eo, i_1c);
model = comsol_fns.create_geometry(model);
model = comsol_fns.add_model_variables(model);
model = comsol_fns.add_materials(model);
model = comsol_fns.create_mesh(model);
model = comsol_fns.create_voltage_probe(model);
model = comsol_fns.add_electrochem_pdes(model, Vo);
model = comsol_fns.add_electrochem_study(model, interval, duration);


% Export geometry png here
model.result.export.create('img1', 'Image');
model.result.export('img1').set('imagetype', 'png');
model.result.export('img1').set('sourcetype', 'geometry');
model.result.export('img1').set('sourceobject', 'geom1');
model.result.export('img1').set('pngfilename', '/Users/SamLy/Desktop/COMSOL API/Matlab/testing.png');
model.result.export('img1').run();
model.result.export().remove('img1');

% % Probably "easier" to entire geometry and recreate it
% model.component('comp1').geom.remove('geom1');

% Try removing the particle sequence instead
model.geom('part1').feature.clear;

% Re-generate ensemble with another porosity
eps2 = 0.65;
% Add particles to existing "Geometry parts"
[circles, model] = comsol_fns.generate_particles(model, ...
    min_r, max_r, clearance, eps2, l_e, h_cell);
porosity = Circle.porosity(circles, l_e, h_cell);
[mean, std] = Circle.particle_stats(circles);

model.component('comp1').geom('geom1').run;

model = comsol_fns.run_electrochem_study(model);

disp('break here');
