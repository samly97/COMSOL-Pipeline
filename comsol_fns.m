classdef comsol_fns
    methods(Static)
        function model = setup_model(path, fname)
            import com.comsol.model.*
            import com.comsol.model.util.*
            
            global PARTICLE_TAG
            
            PARTICLE_TAG = 'part1';
            
            % Create COMSOL file and name it
            model = ModelUtil.create('Model');
            model.modelPath(path);
            model.label(fname);
            
            model.component.create('comp1', true);
            model.component('comp1').geom.create('geom1', 2);
            model.result.table.create('tbl1', 'Table');
            model.component('comp1').mesh.create('mesh1');
            
            % Create "housing" for importing particles
            model.geom.create(PARTICLE_TAG, 'Part', 2);
            model.geom(PARTICLE_TAG).label('Particles');
        end
        
        function [circles, model] = generate_particles(model, ...
                min_r, max_r, clearance, eps, l_e, h_cell)
            % generate_particle fills a rectangular region (l_e, h_cell)
            % with random (uniformly) generated spheres between min_r and
            % max_r. Option to add "clearance" between particles.
            %
            % Inputs:
            % min_r: minimum particle radius (um)
            % max_r: maximum particle radius (um)
            % clearance: spacing between particles (um)
            % eps: desired electrode porosity
            % l_e: length of electrode (um)
            % h_cell: height of electrode (um)
            %
            % Returns
            % circles: array of "Circle" objects, mostly for particle
            % distribution statistics for post-processing.
            % model: the COMSOL model being created.
            
            global PARTICLE_TAG
            
            circles = RSA(min_r, max_r, clearance, eps, l_e, h_cell);
            for i = 1:length(circles)
                circle = circles(i);
                idx = "c" + num2str(i);
                model.geom(PARTICLE_TAG).create(idx, 'Circle');
                model.geom(PARTICLE_TAG).feature(idx).set('pos', [circle.x circle.y]);
                model.geom(PARTICLE_TAG).feature(idx).set('r', circle.R);
            end
            model.geom(PARTICLE_TAG).run;
        end
        
        function model = add_model_parameters(model, h_cell, l_e, l_sep,...
                C_so, C_eo, i_1c)
            % Separator params
            model.param.set('brug', '1.5', 'Bruggeman coefficient');
            model.param.set('eps_sep_l', '0.4', 'Electrolyte Volume fraction in separator');
            
            % Electrode params
            model.param.group.create('default2');
            model.param('default2').set('h_cell', sprintf('%d [um]', h_cell), ...
                'Electrode Height');
            model.param('default2').set('L_pos', sprintf('%d [um]', l_e), ...
                'Positive electrode thickness');
            model.param('default2').set('L_sep', sprintf('%d [um]', l_sep), ...
                'Separator thickness');
            
            % Electrolyte params
            model.param.group.create('par2');
            model.param('par2').set('De', '7.5e-11 * 10 [m^2/s]', 'Salt diffusivity in electrolyte');
            model.param('par2').set('tp', '0.363', 'Lithium ion transference number');
            model.param('par2').set('dlnfdlnc', '0');
            model.param('par2').set('Ce_0', sprintf('%d [mol/m^3]', C_eo), ...
                'Initial salt concentration');
            
            % Particle (NMC) params
            model.param.group.create('par3');
            model.param('par3').set('SOC_max', '0.975', 'Maximum electrode SOC');
            model.param('par3').set('SOC_min', '0', 'Minimum electrode State of Charge');
            model.param('par3').set('Cs_max', '48900 [mol/m^3]', 'Maximum Lithium ion concentration in NMC electrode');
            model.param('par3').set('Ds', '5e-13[m^2/s]', 'Lithium diffusivity in NMC');
            model.param('par3').set('sigma', '100 [S/m]', 'NMC electrical conductivity');
            model.param('par3').set('Cs_0', sprintf('%d [mol/m^3]', C_so), ...
                'Initial Lithium concentration in NMC');
            model.param('par3').set('k', '2e-11 [m/s]', 'Cathodic/Anodic reaction rate constant');
            model.param('par3').set('k_norm', 'k*Cs_max', 'Normalized reaction rate constant');
            model.param('par3').set('alpha', '0.5', 'Charge transfer coefficient');
            model.param('par3').set('sigma_s', '1 [S/m]', 'Porous binder conductivity');
            
            % Operating params
            model.param.group.create('par4');
            model.param('par4').set('T', '20 + 273[K]', 'Cell temperature during operation');
            model.param('par4').set('F', 'F_const', 'Faraday''s constant');
            model.param('par4').set('R', 'R_const', 'Ideal gas constant');
            model.param('par4').set('i_app', sprintf('%f [A/m^2]', i_1c), ...
                'Applied current density');
            model.param('par4').set('LVC', '3.1[V]', 'Lower voltage cutoff');
            
            model.param('default2').label('Geometry Parameters');
            model.param.label('Separator');
            model.param('par2').label('Electrolyte');
            model.param('par3').label('NMC Electrode');
            model.param('par4').label('Operational Parameters');
        end
        
        function model = modify_applied_current(model, i_1c, C)
            % modify_applied_current changes the "i_app" parameter in the
            % "General Parameters" node. This variable is the "top-level"
            % variable, so the experiment will follow this value.
            %
            % model: the COMSOL model object
            % i_1c: the 1C current density for the microstructure (A/m^2)
            % C: discharged C rate (1)
            
            i_app = sprintf('%f [A/m^2]', i_1c * C);
            model.param('par4').set('i_app', i_app);
        end
        
        function model = create_geometry(model)
            % Import particles
            model.component('comp1').geom('geom1').create('pi1', 'PartInstance');
            model.component('comp1').geom('geom1').feature('pi1').set('selkeepnoncontr', false);
            
            % Join particles
            model.component('comp1').geom('geom1').create('uni1', 'Union');
            model.component('comp1').geom('geom1').feature('uni1').label('Particles');
            model.component('comp1').geom('geom1').feature('uni1').set('selresult', true);
            model.component('comp1').geom('geom1').feature('uni1').set('intbnd', false);
            model.component('comp1').geom('geom1').feature('uni1').selection('input').set({'pi1'});
            
            model.component('comp1').geom('geom1').create('r1', 'Rectangle');
            model.component('comp1').geom('geom1').feature('r1').set('pos', [0 0]);
            model.component('comp1').geom('geom1').feature('r1').set('size', {'L_pos' 'h_cell'});
            
            model.component('comp1').geom('geom1').create('par1', 'Partition');
            model.component('comp1').geom('geom1').feature('par1').selection('input').set({'uni1'});
            model.component('comp1').geom('geom1').feature('par1').selection('tool').set({'r1'});
            
            model.component('comp1').geom('geom1').create('boxsel1', 'BoxSelection');
            model.component('comp1').geom('geom1').feature('boxsel1').set('xmin', 0);
            model.component('comp1').geom('geom1').feature('boxsel1').set('xmax', 'L_pos');
            model.component('comp1').geom('geom1').feature('boxsel1').set('ymin', 0);
            model.component('comp1').geom('geom1').feature('boxsel1').set('ymax', 'h_cell');
            model.component('comp1').geom('geom1').feature('boxsel1').set('condition', 'allvertices');
            model.component('comp1').geom('geom1').feature('boxsel1').set('selshow', false);
            
            model.component('comp1').geom('geom1').create('comsel1', 'ComplementSelection');
            model.component('comp1').geom('geom1').feature('comsel1').set('input', {'boxsel1'});
            model.component('comp1').geom('geom1').feature('comsel1').set('selshow', false);
            
            model.component('comp1').geom('geom1').create('del1', 'Delete');
            model.component('comp1').geom('geom1').feature('del1').selection('input').init(2);
            model.component('comp1').geom('geom1').feature('del1').selection('input').named('comsel1');
            
            model.component('comp1').geom('geom1').create('r2', 'Rectangle');
            model.component('comp1').geom('geom1').feature('r2').label('Electrode');
            model.component('comp1').geom('geom1').feature('r2').set('size', {'L_pos' 'h_cell'});
            
            model.component('comp1').geom('geom1').create('dif1', 'Difference');
            model.component('comp1').geom('geom1').feature('dif1').set('keep', true);
            model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'r2'});
            model.component('comp1').geom('geom1').feature('dif1').selection('input2').named('uni1');
            
            model.component('comp1').geom('geom1').create('del2', 'Delete');
            model.component('comp1').geom('geom1').feature('del2').selection('input').init;
            model.component('comp1').geom('geom1').feature('del2').selection('input').set({'r2'});
            
            model.component('comp1').geom('geom1').create('r3', 'Rectangle');
            model.component('comp1').geom('geom1').feature('r3').label('Separator');
            model.component('comp1').geom('geom1').feature('r3').set('selresult', true);
            model.component('comp1').geom('geom1').feature('r3').set('pos', {'-L_sep' '0'});
            model.component('comp1').geom('geom1').feature('r3').set('size', {'L_sep' 'h_cell'});
            
            model.component('comp1').geom('geom1').create('r4', 'Rectangle');
            model.component('comp1').geom('geom1').feature('r4').label('Lithium Foil');
            model.component('comp1').geom('geom1').feature('r4').set('selresult', true);
            model.component('comp1').geom('geom1').feature('r4').set('pos', {'-L_sep - L_pos/10' '0'});
            model.component('comp1').geom('geom1').feature('r4').set('size', {'L_pos/10' 'h_cell'});
            
            % Form union
            model.component('comp1').geom('geom1').run('fin');
            
            model.component('comp1').geom('geom1').create('comsel2', 'ComplementSelection');
            model.component('comp1').geom('geom1').feature('comsel2').label('Porous Electrode Binder');
            model.component('comp1').geom('geom1').feature('comsel2').set('input', {'uni1' 'r3' 'r4'});
            
            model.component('comp1').geom('geom1').create('adjsel1', 'AdjacentSelection');
            model.component('comp1').geom('geom1').feature('adjsel1').label('Particle Boundaries');
            model.component('comp1').geom('geom1').feature('adjsel1').set('input', {'uni1'});
            
            model.component('comp1').geom('geom1').create('adjsel2', 'AdjacentSelection');
            model.component('comp1').geom('geom1').feature('adjsel2').label('Binder Boundaries');
            model.component('comp1').geom('geom1').feature('adjsel2').set('input', {'comsel2'});
            
            model.component('comp1').geom('geom1').create('intsel1', 'IntersectionSelection');
            model.component('comp1').geom('geom1').feature('intsel1').set('entitydim', 1);
            model.component('comp1').geom('geom1').feature('intsel1').label('Particle Surfaces');
            model.component('comp1').geom('geom1').feature('intsel1').set('input', {'adjsel1' 'adjsel2'});
            
            model.component('comp1').geom('geom1').create('unisel1', 'UnionSelection');
            model.component('comp1').geom('geom1').feature('unisel1').label('Separator + Binder');
            model.component('comp1').geom('geom1').feature('unisel1').set('input', {'r3' 'comsel2'});
            
            model.component('comp1').geom('geom1').create('unisel2', 'UnionSelection');
            model.component('comp1').geom('geom1').feature('unisel2').label('Foil + Particles');
            model.component('comp1').geom('geom1').feature('unisel2').set('input', {'uni1' 'r4'});
            
            model.component('comp1').geom('geom1').create('unisel3', 'UnionSelection');
            model.component('comp1').geom('geom1').feature('unisel3').label('Foil + Particles + Binder');
            model.component('comp1').geom('geom1').feature('unisel3').set('input', {'uni1' 'r4' 'comsel2'});
            
            model.component('comp1').geom('geom1').create('boxsel3', 'BoxSelection');
            model.component('comp1').geom('geom1').feature('boxsel3').set('entitydim', 1);
            model.component('comp1').geom('geom1').feature('boxsel3').label('Current Collector');
            model.component('comp1').geom('geom1').feature('boxsel3').set('xmin', 'L_pos');
            model.component('comp1').geom('geom1').feature('boxsel3').set('xmax', 'L_pos');
            model.component('comp1').geom('geom1').feature('boxsel3').set('condition', 'allvertices');
            
            model.component('comp1').geom('geom1').create('boxsel2', 'BoxSelection');
            model.component('comp1').geom('geom1').feature('boxsel2').label('All Domains');
            model.component('comp1').geom('geom1').feature('boxsel2').set('condition', 'allvertices');
            model.component('comp1').geom('geom1').feature('boxsel2').set('selshow', false);
            
            model.component('comp1').geom('geom1').create('adjsel3', 'AdjacentSelection');
            model.component('comp1').geom('geom1').feature('adjsel3').label('All External Boundaries');
            model.component('comp1').geom('geom1').feature('adjsel3').set('input', {'boxsel2'});
            
            % Finalize geometry and run
            model.component('comp1').geom('geom1').run;
        end
        
        function model = add_model_variables(model)
            model.component('comp1').variable.create('var1');
            model.component('comp1').variable('var1').set('De_eff', 'eps_sep_l^brug*De');
            model.component('comp1').variable('var1').set('kap_eff', 'eps_sep_l^brug*mat2.ionc.sigmal_int1(C_e/Ce_0)*1[S/m]');
            model.component('comp1').variable('var1').set('kap_eff_D', 'kap_eff*2*R*T/F*(1+dlnfdlnc)*(1-tp)');
            model.component('comp1').variable('var1').selection.named('geom1_r3_dom');
            model.component('comp1').variable.create('var2');
            model.component('comp1').variable('var2').set('kap', 'mat2.ionc.sigmal_int1(C_e/Ce_0)*1[S/m]');
            model.component('comp1').variable('var2').set('kap_D', 'kap*2*R*T/F*(1+dlnfdlnc)*(1-tp)');
            model.component('comp1').variable('var2').selection.named('geom1_comsel2');
            model.component('comp1').variable.create('var3');
            model.component('comp1').variable('var3').set('soc', 'C_s/Cs_max');
            model.component('comp1').variable('var3').set('eta', 'phi_s - phi_l - mat1.elpot.Eeq_int1(soc)');
            model.component('comp1').variable('var3').set('j0', 'k_norm*((Cs_max - C_s)/Cs_max)^alpha*(C_s/Cs_max)^alpha*(C_e/1[mol/m^3])^alpha');
            model.component('comp1').variable('var3').set('j_bv', 'j0*(exp(alpha*F*eta/(R*T)) - exp(-alpha*F*eta/(R*T)))', 'Butler-Volmer reaction rate in molar form (mol/m^2/s)');
            
            model.component('comp1').variable('var1').label('Separator Variables');
            model.component('comp1').variable('var2').label('Electrolyte Variables');
            model.component('comp1').variable('var3').label('Reaction Variables');
        end
        
        function model = add_materials(model)
        % Material 1 - NMC Particles
        % Material 2 - LiPF6 in 1:2 EC:DMC and p(VdF-HFP)
        model = comsol_fns.h_make_material_shell(model);
        model = comsol_fns.h_add_material_nmc(model);
        model = comsol_fns.h_add_material_electrolyte(model);
        end
        
        function model = h_make_material_shell(model)
            model.component('comp1').material.create('mat1', 'Common');
            model.component('comp1').material('mat1').propertyGroup.create('ElectrodePotential', 'Equilibrium potential');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').func.create('int1', 'Interpolation');
            model.component('comp1').material('mat1').propertyGroup.create('OperationalSOC', 'Operational electrode state-of-charge');
            
            model.component('comp1').material.create('mat2', 'Common');
            model.component('comp1').material('mat2').propertyGroup.create('ElectrolyteConductivity', 'Electrolyte conductivity');
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteConductivity').func.create('int1', 'Interpolation');
            model.component('comp1').material('mat2').propertyGroup.create('SpeciesProperties', 'Species properties');
            model.component('comp1').material('mat2').propertyGroup.create('ElectrolyteSaltConcentration', 'Electrolyte salt concentration');
        end
        
        function model = h_add_material_nmc(model)
            model.component('comp1').material('mat1').label('NMC 111 Electrode, LiNi0.33Mn0.33Co0.33O2 (Positive, Li-ion Battery)');
            model.component('comp1').material('mat1').propertyGroup('def').set('diffusion', {'5e-13[m^2/s]' '0' '0' '0' '5e-13[m^2/s]' '0' '0' '0' '5e-13[m^2/s]'});
            model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:diffusion', 'W. Zheng, M. Shui, J. Shu, S. Gao, D. Xu, L. Chen, L. Feng and Y. Ren, " GITT studies on oxide cathode LiNi1/3Co1/3Mn1/3O2 synthesized by citric acid assisted high-energy ball milling", Bull. Mater. Sci., vol. 36, p. A495, 2013');
            model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'100[S/m]' '0' '0' '0' '100[S/m]' '0' '0' '0' '100[S/m]'});
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').func('int1').set('funcname', 'Eeq_int1');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').func('int1').set('table', {'0' '4.44';  ...
                '0.032' '4.34';  ...
                '0.102' '4.23';  ...
                '0.187' '4.13';  ...
                '0.289' '4.025';  ...
                '0.38' '3.945';  ...
                '0.543' '3.835';  ...
                '0.775' '3.71';  ...
                '0.872' '3.62';  ...
                '0.925' '3.51';  ...
                '0.943' '3.42';  ...
                '0.957' '3.30';  ...
                '0.966' '3.165';  ...
                '0.970' '3.02';  ...
                '0.972' '2.90';  ...
                '0.975' '2.688'});
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').func('int1').set('interp', 'piecewisecubic');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').func('int1').set('extrap', 'linear');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').func('int1').set('fununit', 'V');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('Eeq', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('dEeqdT', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('cEeqref', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('Eeq', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('dEeqdT', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('cEeqref', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('Eeq', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('dEeqdT', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('cEeqref', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('Eeq', 'Eeq_int1(soc)+dEeqdT*(T-298[K])');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('INFO_PREFIX:Eeq', 'W. Zheng, M. Shui, J. Shu, S. Gao, D. Xu, L. Chen, L. Feng and Y. Ren, " GITT studies on oxide cathode LiNi1/3Co1/3Mn1/3O2 synthesized by citric acid assisted high-energy ball milling", Bull. Mater. Sci., vol. 36, p. A495, 2013');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('dEeqdT', '-10[J/mol/K]/F_const');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('INFO_PREFIX:dEeqdT', 'V Viswanathan, D Choi, D Wang, W Xu, S Towne, R Williford, JG Zhang, J Liu and Z Yang "Effect of entropy change on lithium intercalation in cathodes and anodes on Li-ion battery thermal management", Journal of Power Sources 195 (2010) 3720-3729');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('cEeqref', '49000[mol/m^3]');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('INFO_PREFIX:cEeqref', 'W. Zheng, M. Shui, J. Shu, S. Gao, D. Xu, L. Chen, L. Feng and Y. Ren, " GITT studies on oxide cathode LiNi1/3Co1/3Mn1/3O2 synthesized by citric acid assisted high-energy ball milling", Bull. Mater. Sci., vol. 36, p. A495, 2013');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').set('soc', 'c/cEeqref');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').descr('soc', '');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').addInput('concentration');
            model.component('comp1').material('mat1').propertyGroup('ElectrodePotential').addInput('temperature');
            model.component('comp1').material('mat1').propertyGroup('OperationalSOC').set('socmax', '');
            model.component('comp1').material('mat1').propertyGroup('OperationalSOC').set('socmin', '');
            model.component('comp1').material('mat1').propertyGroup('OperationalSOC').set('socmax', '');
            model.component('comp1').material('mat1').propertyGroup('OperationalSOC').set('socmin', '');
            model.component('comp1').material('mat1').propertyGroup('OperationalSOC').set('socmax', '0.975');
            model.component('comp1').material('mat1').propertyGroup('OperationalSOC').set('socmin', '0');
        end
        
        function model = h_add_material_electrolyte(model)
            model.component('comp1').material('mat2').label('LiPF6 in 1:2 EC:DMC and p(VdF-HFP) (Polymer electrolyte, Li-ion Battery)');
            model.component('comp1').material('mat2').comments(['\n']);
            model.component('comp1').material('mat2').propertyGroup('def').set('diffusion', {'7.5e-11[m^2/s]' '0' '0' '0' '7.5e-11[m^2/s]' '0' '0' '0' '7.5e-11[m^2/s]'});
            model.component('comp1').material('mat2').propertyGroup('def').set('INFO_PREFIX:diffusion', ['M. Doyle, J. Newman, A.S. Gozdz, C.N. Schmutz, and J.M. Tarascon, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Comparison of Modeling Predictions with Experimental Data from Plastic Lithium Ion Cells,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Electrochem. Soc., vol. 143, p. 1890, 1996.\n']);
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteConductivity').func('int1').set('funcname', 'sigmal_int1');
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteConductivity').func('int1').set('table', {'0' '0.0108';  ...
                '0.2000' '0.1259';  ...
                '0.4000' '0.2055';  ...
                '0.6000' '0.2553';  ...
                '0.8000' '0.2810';  ...
                '1.0000' '0.2873';  ...
                '1.2000' '0.2788';  ...
                '1.4000' '0.2595';  ...
                '1.6000' '0.2331';  ...
                '1.8000' '0.2027';  ...
                '2.0000' '0.1710';  ...
                '2.200' '0.1403';  ...
                '2.4000' '0.1123';  ...
                '2.6000' '0.0885';  ...
                '2.8000' '0.0696';  ...
                '3.0000' '0.0563'});
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteConductivity').set('sigmal', {'sigmal_int1(c/c_ref)' '0' '0' '0' 'sigmal_int1(c/c_ref)' '0' '0' '0' 'sigmal_int1(c/c_ref)'});
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteConductivity').set('INFO_PREFIX:sigmal', ['M. Doyle, J. Newman, A.S. Gozdz, C.N. Schmutz, and J.M. Tarascon, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Comparison of Modeling Predictions with Experimental Data from Plastic Lithium Ion Cells,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Electrochem. Soc., vol. 143, p. 1890, 1996.\n']);
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteConductivity').set('c_ref', '1000[mol/m^3]');
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteConductivity').descr('c_ref', '');
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteConductivity').addInput('concentration');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('transpNum', '');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('fcl', '');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('transpNum', '');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('fcl', '');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('transpNum', '');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('fcl', '');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('transpNum', '0.363');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('INFO_PREFIX:transpNum', ['M. Doyle, J. Newman, A.S. Gozdz, C.N. Schmutz, and J.M. Tarascon, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Comparison of Modeling Predictions with Experimental Data from Plastic Lithium Ion Cells,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Electrochem. Soc., vol. 143, p. 1890, 1996.\n']);
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('fcl', '0');
            model.component('comp1').material('mat2').propertyGroup('SpeciesProperties').set('INFO_PREFIX:fcl', ['M. Doyle, J. Newman, A.S. Gozdz, C.N. Schmutz, and J.M. Tarascon, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Comparison of Modeling Predictions with Experimental Data from Plastic Lithium Ion Cells,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Electrochem. Soc., vol. 143, p. 1890, 1996.\n']);
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteSaltConcentration').identifier('cElsalt');
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteSaltConcentration').set('cElsalt', '1000[mol/m^3]');
            model.component('comp1').material('mat2').propertyGroup('ElectrolyteSaltConcentration').set('INFO_PREFIX:cElsalt', ['M. Doyle, J. Newman, A.S. Gozdz, C.N. Schmutz, and J.M. Tarascon, ' native2unicode(hex2dec({'20' '1c'}), 'unicode') 'Comparison of Modeling Predictions with Experimental Data from Plastic Lithium Ion Cells,' native2unicode(hex2dec({'20' '1d'}), 'unicode') ' J. Electrochem. Soc., vol. 143, p. 1890, 1996.\n']);
        end
        
        function model = create_mesh(model)
            model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
            model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
            model.component('comp1').mesh('mesh1').create('bl1', 'BndLayer');
            model.component('comp1').mesh('mesh1').feature('ftri2').selection.named('geom1_comsel2');
            model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
            model.component('comp1').mesh('mesh1').feature('bl1').selection.named('geom1_uni1_dom');
            model.component('comp1').mesh('mesh1').feature('bl1').create('blp', 'BndLayerProp');
            model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').selection.named('geom1_intsel1');
            
            model.component('comp1').mesh('mesh1').feature('size').set('hauto', 6);
            model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 9);
            model.component('comp1').mesh('mesh1').feature('bl1').feature('blp').set('blnlayers', 3);
            model.component('comp1').mesh('mesh1').run;
        end
        
        function model = create_voltage_probe(model)
        % Boundary probe to measure cell voltage. Can use the
        % Lower-Cutoff-Voltage to terminal simulation.
        model.component('comp1').probe.create('bnd1', 'Boundary');
        model.component('comp1').probe('bnd1').selection.named('geom1_boxsel3');
        
        % Create probe results table
        model.result.table('tbl1').label('Probe Table 1');
        
        model.component('comp1').probe('bnd1').set('expr', 'phi_s');
        model.component('comp1').probe('bnd1').set('unit', 'V');
        model.component('comp1').probe('bnd1').set('descr', 'Dependent variable phi_s');
        model.component('comp1').probe('bnd1').set('table', 'tbl1');
        model.component('comp1').probe('bnd1').set('window', 'window5');
        end
        
        function model = add_electrochem_pdes(model, Vo)
            model = comsol_fns.h_electrolyte_conc_interface(model);
            model = comsol_fns.h_electrolyte_conc_eqns(model);
            model = comsol_fns.h_solid_conc_interface(model);
            model = comsol_fns.h_solid_conc_eqns(model);
            model = comsol_fns.h_ionic_potential_interface(model);
            model = comsol_fns.h_ionic_potential_eqns(model);
            model = comsol_fns.h_solid_potential_interface(model);
            model = comsol_fns.h_solid_potential_eqns(model, Vo);
            model = comsol_fns.h_discharge_event(model); 
        end
        
        function model = h_electrolyte_conc_interface(model)
            model.component('comp1').physics.create('C_e', 'GeneralFormPDE', 'geom1');
            model.component('comp1').physics('C_e').identifier('C_e');
            model.component('comp1').physics('C_e').field('dimensionless').field('C_e');
            model.component('comp1').physics('C_e').field('dimensionless').component({'C_e'});
            model.component('comp1').physics('C_e').prop('Units').set('DependentVariableQuantity', 'concentration');
            
            % Select separator and binder as domain
            model.component('comp1').physics('C_e').selection.named('geom1_unisel1');
            
            % Select binder as domain in a separate PDE
            model.component('comp1').physics('C_e').create('gfeq2', 'GeneralFormPDE', 2);
            model.component('comp1').physics('C_e').feature('gfeq2').selection.named('geom1_comsel2');
            
            % Set up boundary conditions
            model.component('comp1').physics('C_e').create('flux1', 'FluxBoundary', 1);
            model.component('comp1').physics('C_e').feature('flux1').selection.named('geom1_intsel1');
            model.component('comp1').physics('C_e').create('flux2', 'FluxBoundary', 1);
            model.component('comp1').physics('C_e').feature('flux2').selection.set([4]);
        end
        
        function model = h_electrolyte_conc_eqns(model)
            model.component('comp1').physics('C_e').label('Electrolyte Diffusion');
            model.component('comp1').physics('C_e').prop('Units').set('SourceTermQuantity', 'reactionrate');
            
            % Initial value for concentration
            model.component('comp1').physics('C_e').feature('init1').set('C_e', 'Ce_0');
            
            % Separator PDE Config
            model.component('comp1').physics('C_e').feature('gfeq1').set('f', 0);
            model.component('comp1').physics('C_e').feature('gfeq1').set('Ga', {'-De_eff*C_ex + tp/F*(-kap_eff*phi_lx + kap_eff_D*1/C_e*C_ex)' '-De_eff*C_ey + tp/F*(-kap_eff*phi_ly + kap_eff_D*1/C_e*C_ey)'});
            model.component('comp1').physics('C_e').feature('gfeq1').set('da', 'eps_sep_l');
            model.component('comp1').physics('C_e').feature('gfeq1').label('Separator');
            
            % Porous Binder PDE config
            model.component('comp1').physics('C_e').feature('gfeq2').set('f', 0);
            model.component('comp1').physics('C_e').feature('gfeq2').set('Ga', {'-De*C_ex + tp/F*(-kap*phi_lx + kap_D*1/C_e*C_ex)' '-De*C_ey + tp/F*(-kap*phi_ly + kap_D*1/C_e*C_ey)'});
            model.component('comp1').physics('C_e').feature('gfeq2').label('Electrolyte');
        
            % Reaction on particle surfaces
            model.component('comp1').physics('C_e').feature('flux1').set('g', 'j_bv');
            model.component('comp1').physics('C_e').feature('flux1').label('Flux/Source - Electrode Reaction');

            % Flux at lithium foil
            model.component('comp1').physics('C_e').feature('flux2').set('g', 'i_app/F');
            model.component('comp1').physics('C_e').feature('flux2').label('Flux/Source - Applied Current');
        end
        
        function model = h_solid_conc_interface(model)
            model.component('comp1').physics.create('C_s', 'GeneralFormPDE', 'geom1');
            model.component('comp1').physics('C_s').identifier('C_s');
            model.component('comp1').physics('C_s').field('dimensionless').field('C_s');
            model.component('comp1').physics('C_s').field('dimensionless').component({'C_s'});
            model.component('comp1').physics('C_s').prop('Units').set('DependentVariableQuantity', 'concentration');
        
            % Select NMC as domain
            model.component('comp1').physics('C_s').selection.named('geom1_uni1_dom');
            
            model.component('comp1').physics('C_s').create('flux1', 'FluxBoundary', 1);
            model.component('comp1').physics('C_s').feature('flux1').selection.named('geom1_intsel1');
        end
        
        function model = h_solid_conc_eqns(model)
            model.component('comp1').physics('C_s').label('Solid Diffusion');
            model.component('comp1').physics('C_s').prop('Units').set('SourceTermQuantity', 'reactionrate');
            
            % Solid diffusion equation
            model.component('comp1').physics('C_s').feature('gfeq1').set('f', 0);
            model.component('comp1').physics('C_s').feature('gfeq1').set('Ga', {'-Ds*C_sx' '-Ds*C_sy'});
            
            % initial conc and Rxn on NMC surface
            model.component('comp1').physics('C_s').feature('init1').set('C_s', 'Cs_0');
            model.component('comp1').physics('C_s').feature('flux1').set('g', '-j_bv');
        end
        
        function model = h_ionic_potential_interface(model)
            model.component('comp1').physics.create('phi_l', 'GeneralFormPDE', 'geom1');
            model.component('comp1').physics('phi_l').identifier('phi_l');
            model.component('comp1').physics('phi_l').field('dimensionless').field('phi_l');
            model.component('comp1').physics('phi_l').field('dimensionless').component({'phi_l'});
            model.component('comp1').physics('phi_l').prop('Units').set('DependentVariableQuantity', 'thermoneutral_voltage');
            
            % Porous binder and separator as domain
            model.component('comp1').physics('phi_l').selection.named('geom1_unisel1');
            
            % Define PDE on porous binder region
            model.component('comp1').physics('phi_l').create('gfeq2', 'GeneralFormPDE', 2);
            model.component('comp1').physics('phi_l').feature('gfeq2').selection.named('geom1_comsel2');
            
            % Flux at lithium foil
            model.component('comp1').physics('phi_l').create('flux1', 'FluxBoundary', 1);
            model.component('comp1').physics('phi_l').feature('flux1').selection.set([4]);
            
            % Flux at NMC particle surfaces
            model.component('comp1').physics('phi_l').create('flux2', 'FluxBoundary', 1);
            model.component('comp1').physics('phi_l').feature('flux2').selection.named('geom1_intsel1');

            % Ground node at Lithium foil
            model.component('comp1').physics('phi_l').create('cons1', 'Constraint', 1);
            model.component('comp1').physics('phi_l').feature('cons1').selection.set([4]);
        end
        
        function model = h_ionic_potential_eqns(model)
            model.component('comp1').physics('phi_l').label('Electrolyte Potential');
            model.component('comp1').physics('phi_l').prop('Units').set('SourceTermQuantity', 'currentsource');
            
            % Separator PDE
            model.component('comp1').physics('phi_l').feature('gfeq1').set('f', 0);
            model.component('comp1').physics('phi_l').feature('gfeq1').set('Ga', {'-kap_eff*phi_lx + kap_eff_D*1/C_e*C_ex' '-kap_eff*phi_ly + kap_eff_D*1/C_e*C_ey'});
            model.component('comp1').physics('phi_l').feature('gfeq1').set('da', 0);
            model.component('comp1').physics('phi_l').feature('gfeq1').label('Separator');
            
            model.component('comp1').physics('phi_l').feature('init1').label('Initial Values ');
            
            % Electrolyte PDE
            model.component('comp1').physics('phi_l').feature('gfeq2').set('f', 0);
            model.component('comp1').physics('phi_l').feature('gfeq2').set('Ga', {'-kap*phi_lx + kap_D*1/C_e*C_ex' '-kap*phi_ly + kap_D*1/C_e*C_ey'});
            model.component('comp1').physics('phi_l').feature('gfeq2').set('da', 0);
            model.component('comp1').physics('phi_l').feature('gfeq2').label('Electrolyte');
            
            % Current at Lithium foil
            model.component('comp1').physics('phi_l').feature('flux1').set('g', 'i_app');
            model.component('comp1').physics('phi_l').feature('flux1').label('Flux/Source - Applied Current');
            
            % Current at NMC particle surfaces
            model.component('comp1').physics('phi_l').feature('flux2').set('g', 'j_bv*F');
            model.component('comp1').physics('phi_l').feature('flux2').label('Flux/Source - Reaction');
            
            % Ground node at Lithium foil
            model.component('comp1').physics('phi_l').feature('cons1').set('R', 'phi_l');
            model.component('comp1').physics('phi_l').feature('cons1').label('Ground Node');
        end
        
        function model = h_solid_potential_interface(model)
            model.component('comp1').physics.create('phi_s', 'GeneralFormPDE', 'geom1');
            model.component('comp1').physics('phi_s').identifier('phi_s');
            model.component('comp1').physics('phi_s').field('dimensionless').field('phi_s');
            model.component('comp1').physics('phi_s').field('dimensionless').component({'phi_s'});
            model.component('comp1').physics('phi_s').prop('Units').set('DependentVariableQuantity', 'thermoneutral_voltage');
            
            % Select Lithium foil, particles, and binder
            model.component('comp1').physics('phi_s').selection.named('geom1_unisel3');
            
            % Initial value - NMC Particles
            model.component('comp1').physics('phi_s').create('init3', 'init', 2);
            model.component('comp1').physics('phi_s').feature('init3').selection.named('geom1_uni1_dom');
            
            % Initial value - Lithium foil
            model.component('comp1').physics('phi_s').create('init2', 'init', 2);
            model.component('comp1').physics('phi_s').feature('init2').selection.named('geom1_r4_dom');
            
            % Lithium foil PDE
            model.component('comp1').physics('phi_s').create('gfeq2', 'GeneralFormPDE', 2);
            model.component('comp1').physics('phi_s').feature('gfeq2').selection.named('geom1_r4_dom');
            
            % Porous binder PDE
            model.component('comp1').physics('phi_s').create('gfeq3', 'GeneralFormPDE', 2);
            model.component('comp1').physics('phi_s').feature('gfeq3').selection.named('geom1_comsel2');
            
            % Current at cathode current collector
            model.component('comp1').physics('phi_s').create('flux1', 'FluxBoundary', 1);
            model.component('comp1').physics('phi_s').feature('flux1').selection.named('geom1_boxsel3');
            
            % Current at NMC particle surfaces
            model.component('comp1').physics('phi_s').create('flux2', 'FluxBoundary', 1);
            model.component('comp1').physics('phi_s').feature('flux2').selection.named('geom1_intsel1');
            
            % Ground node at Lithium foil
            model.component('comp1').physics('phi_s').create('cons1', 'Constraint', 1);
            model.component('comp1').physics('phi_s').feature('cons1').selection.set([4]);
        end
        
        function model = h_solid_potential_eqns(model, Vo)
            model.component('comp1').physics('phi_s').label('Solid Potential');
            model.component('comp1').physics('phi_s').prop('Units').set('SourceTermQuantity', 'currentsource');
            
            % NMC Particles PDE
            model.component('comp1').physics('phi_s').feature('gfeq1').set('f', 0);
            model.component('comp1').physics('phi_s').feature('gfeq1').set('Ga', {'-sigma*phi_sx' '-sigma*phi_sy'});
            model.component('comp1').physics('phi_s').feature('gfeq1').set('da', 0);
            model.component('comp1').physics('phi_s').feature('gfeq1').label('NMC Particles');
            
            % Initial voltage - Porous Binder
            model.component('comp1').physics('phi_s').feature('init1').set('phi_s', Vo);
            model.component('comp1').physics('phi_s').feature('init1').label('Initial Values - Binder');
            
            % Initial voltage - NMC Particles
            model.component('comp1').physics('phi_s').feature('init3').set('phi_s', 4.29);
            model.component('comp1').physics('phi_s').feature('init3').label('Initial Values - NMC');
            
            % Lithium Foil PDE
            model.component('comp1').physics('phi_s').feature('gfeq2').set('f', 0);
            model.component('comp1').physics('phi_s').feature('gfeq2').set('Ga', {'-sigma*phi_sx' '-sigma*phi_sy'});
            model.component('comp1').physics('phi_s').feature('gfeq2').set('da', 0);
            model.component('comp1').physics('phi_s').feature('gfeq2').label('Lithium Foil PDE');
            
            % Porous Binder PDE
            model.component('comp1').physics('phi_s').feature('gfeq3').set('f', 0);
            model.component('comp1').physics('phi_s').feature('gfeq3').set('Ga', {'-sigma_s*phi_sx' '-sigma_s*phi_sy'});
            model.component('comp1').physics('phi_s').feature('gfeq3').set('da', 0);
            model.component('comp1').physics('phi_s').feature('gfeq3').label('Binder PDE');
            
            model.component('comp1').physics('phi_s').feature('init2').label('Initial Values - Foil');
            
            % Applied current at NMC current collector
            model.component('comp1').physics('phi_s').feature('flux1').set('g', '-i_app');
            model.component('comp1').physics('phi_s').feature('flux1').label('Flux/Source - Applied Current');
            
            % Reaction at NMC particle surfaces
            model.component('comp1').physics('phi_s').feature('flux2').set('g', '-j_bv*F');
            model.component('comp1').physics('phi_s').feature('flux2').label('Flux/Source - Reaction');
            
            % Ground node at Lithium Foil
            model.component('comp1').physics('phi_s').feature('cons1').set('R', 'phi_s');
            model.component('comp1').physics('phi_s').feature('cons1').label('Electrical Ground');
        end
        
        function model = h_discharge_event(model) 
            model.component('comp1').physics.create('ev', 'Events', 'geom1');
            model.component('comp1').physics('ev').create('impl1', 'ImplicitEvent', -1);
            model.component('comp1').physics('ev').create('is1', 'IndicatorStates', -1);
            
            
            model.component('comp1').physics('ev').feature('impl1').set('condition', 'STOP_DCH < 0');
            model.component('comp1').physics('ev').feature('is1').set('indDim', 'STOP_DCH');
            model.component('comp1').physics('ev').feature('is1').set('g', 'bnd1 - LVC');
            model.component('comp1').physics('ev').feature('is1').set('dimInit', 0);
        end
        
        function model = add_electrochem_study(model) 
        % Adds a stationary and time-dependent study to the COMSOL model to
        % simulate heterogeneous study of Lithium-ion battery operation in
        % 2D.
            model = comsol_fns.h_add_studies(model);
            model = comsol_fns.h_def_solver_settings(model);
            model = comsol_fns.h_probe_results_settings(model);
        end
        
        function model = setup_electrochem_study_duration(model, step, stop)
            % Specifies the interval of data collection and stop time for
            % the electrochemical study.
            
            % step: time-step intervals to get data for (s)
            % stop: how long to run the study for (s)
            
            % Study discharge from 0 - 'stop' with 'step' as the interval
            % of data collection
            discharge_time_data = sprintf('range(0,%d,%d)', step, stop);
            model.study('std1').feature('time').set('tlist', ...
                discharge_time_data);
            
            model = comsol_fns.h_more_solver_settings(model, ...
                discharge_time_data);
        end
        
        function model = h_add_studies(model)
            model.study.create('std1');
            model.study('std1').create('stat', 'Stationary');
            model.study('std1').create('time', 'Transient');
            model.study('std1').feature('stat').set('activate', {'C_e' 'off' ...
                'C_s' 'off' 'phi_l' 'on' 'phi_s' 'on' 'ev' 'off'  ...
                'tds' 'off' 'frame:spatial1' 'on' 'frame:material1' 'on'});
            model.study('std1').feature('time').activate('tds', false);
        end
        
        function model = h_def_solver_settings(model)
            model.sol.create('sol1');
            model.sol('sol1').study('std1');
            model.sol('sol1').attach('std1');
            model.sol('sol1').create('st1', 'StudyStep');
            model.sol('sol1').create('v1', 'Variables');
            model.sol('sol1').create('s1', 'Stationary');
            model.sol('sol1').create('su1', 'StoreSolution');
            model.sol('sol1').create('st2', 'StudyStep');
            model.sol('sol1').create('v2', 'Variables');
            model.sol('sol1').create('t1', 'Time');
            model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
            model.sol('sol1').feature('s1').feature.remove('fcDef');
            model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
            model.sol('sol1').feature('t1').create('se1', 'Segregated');
            model.sol('sol1').feature('t1').create('st1', 'StopCondition');
            model.sol('sol1').feature('t1').feature('se1').create('ev1', 'SegregatedStep');
            model.sol('sol1').feature('t1').feature('se1').create('C_s1', 'SegregatedStep');
            model.sol('sol1').feature('t1').feature('se1').create('C_e1', 'SegregatedStep');
            model.sol('sol1').feature('t1').feature.remove('fcDef');
        end
        
        function model = h_probe_results_settings(model)
            model.result.dataset.create('avh1', 'Average');
            model.result.dataset('dset2').set('probetag', 'bnd1');
            model.result.dataset('dset2').set('solution', 'sol1');
            model.result.dataset('avh1').set('probetag', 'bnd1');
            model.result.dataset('avh1').set('data', 'dset2');
            model.result.dataset('avh1').selection.geom('geom1', 1);
            model.result.dataset('avh1').selection.named('geom1_boxsel3');
            model.result.numerical.create('pev1', 'EvalPoint');
            model.result.numerical('pev1').set('probetag', 'bnd1');
            model.result.create('pg5', 'PlotGroup1D');
            model.result('pg5').set('probetag', 'window5_default');
            model.result('pg5').create('tblp1', 'Table');
            model.result('pg5').feature('tblp1').set('probetag', 'bnd1');
            
            model.component('comp1').probe('bnd1').genResult([]);
            
            model.result('pg1').tag('pg5');
        end
        
        function model = h_more_solver_settings(model, study_duration)
            model.sol('sol1').attach('std1');
            model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
            model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
            model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
            model.sol('sol1').feature('s1').set('stol', '1E-4');
            model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
            model.sol('sol1').feature('s1').feature('dDef').set('linsolver', 'pardiso');
            model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
            model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
            model.sol('sol1').feature('su1').label('Solution Store 1.1');
            model.sol('sol1').feature('st2').label('Compile Equations: Time Dependent');
            model.sol('sol1').feature('st2').set('studystep', 'time');
            model.sol('sol1').feature('v2').label('Dependent Variables 2.1');
            model.sol('sol1').feature('v2').set('control', 'user');
            model.sol('sol1').feature('v2').set('initmethod', 'sol');
            model.sol('sol1').feature('v2').set('initsol', 'sol1');
            model.sol('sol1').feature('v2').set('initsoluse', 'sol2');
            model.sol('sol1').feature('v2').set('solnum', 'auto');
            model.sol('sol1').feature('v2').set('resscalemethod', 'manual');
            model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
            model.sol('sol1').feature('v2').set('notsol', 'sol1');
            model.sol('sol1').feature('v2').set('notsolnum', 'auto');
            
            % Not sure what this line does, but let's replace it with the
            % "generic" duration to be safe
            model.sol('sol1').feature('v2').set('clist', {study_duration '0.001[s]'});
            
            model.sol('sol1').feature('t1').label('Time-Dependent Solver 1.1');
            model.sol('sol1').feature('t1').set('control', 'time');
            
            % This line here, duration
            model.sol('sol1').feature('t1').set('tlist', study_duration);
            
            model.sol('sol1').feature('t1').set('initialstepbdfactive', true);
            model.sol('sol1').feature('t1').set('stabcntrl', true);
            model.sol('sol1').feature('t1').feature('dDef').label('Direct 1');
            model.sol('sol1').feature('t1').feature('dDef').set('linsolver', 'pardiso');
            model.sol('sol1').feature('t1').feature('aDef').label('Advanced 1');
            model.sol('sol1').feature('t1').feature('fc1').label('Fully Coupled 1.1');
            model.sol('sol1').feature('t1').feature('se1').label('Segregated 1.1');
            model.sol('sol1').feature('t1').feature('se1').set('segstabacc', 'segaacc');
            model.sol('sol1').feature('t1').feature('se1').feature('ssDef').label('Segregated Step 2');
            model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('segvar', {'comp1_phi_l' 'comp1_phi_s'});
            model.sol('sol1').feature('t1').feature('se1').feature('ssDef').set('subjtech', 'onevery');
            model.sol('sol1').feature('t1').feature('se1').feature('ev1').label('Segregated Step 1.1');
            model.sol('sol1').feature('t1').feature('se1').feature('ev1').set('segvar', {'comp1_ev_is1_indDim'});
            model.sol('sol1').feature('t1').feature('se1').feature('C_s1').label('Segregated Step 1a 1');
            model.sol('sol1').feature('t1').feature('se1').feature('C_s1').set('segvar', {'comp1_C_s'});
            model.sol('sol1').feature('t1').feature('se1').feature('C_e1').label('Segregated Step 1b 1');
            model.sol('sol1').feature('t1').feature('se1').feature('C_e1').set('segvar', {'comp1_C_e'});
            model.sol('sol1').feature('t1').feature('st1').label('Stop Condition 1.1');
            model.sol('sol1').feature('t1').feature('st1').set('eventstopActive', {'on'});
            model.sol('sol1').feature('t1').feature('st1').set('stopcondwarn', false);
        end
        
        function model = run_electrochem_study(model)
            model.sol('sol1').runAll;
        end
        
        function model = export_geometry_pic(model, filepath)
        % export_geometry_pic exports the picture of the cell and then
        % removes the image from the "Results > Export" nodes in COMSOL.
        %
        % model: comsol.model object
        % filepath: where to save the exported image
        model.result.export.create('img1', 'Image');
        model.result.export('img1').set('zoomextents', 'on');
        model.result.export('img1').set('imagetype', 'png');
        model.result.export('img1').set('sourcetype', 'geometry');
        model.result.export('img1').set('sourceobject', 'geom1');
        model.result.export('img1').set('pngfilename', filepath);
        model.result.export('img1').run();
        model.result.export().remove('img1');
        end
        
        function model = add_dilute_transport(model)
            % h_add_dilute_transport adds the "Transport of Dilute Species"
            % interface to "component 1". Changes the diffusion coefficient
            % to 1, adds the initial and boundary conditions.
            model.component('comp1').physics.create('tds', 'DilutedSpecies', {'c'});
            model.component('comp1').physics('tds').selection.named('geom1_comsel2');
            model.component('comp1').physics('tds').label('Tortuosity Measurement');
            model.component('comp1').physics('tds').prop('TransportMechanism').set('Convection', false);
            
            model.component('comp1').physics('tds').feature('cdm1').set('D_c', {'1[m^2/s]' '0' '0' '0' '1[m^2/s]' '0' '0' '0' '1[m^2/s]'});
            model.component('comp1').physics('tds').create('conc1', 'Concentration', 1);
            model.component('comp1').physics('tds').create('conc2', 'Concentration', 1);
            model.component('comp1').physics('tds').feature('conc1').label('Concentration - Inlet');
            model.component('comp1').physics('tds').feature('conc1').selection.set([7]);
            model.component('comp1').physics('tds').feature('conc2').label('Concentration - Outlet');
            model.component('comp1').physics('tds').feature('conc2').selection.named('geom1_boxsel3');
            model.component('comp1').physics('tds').feature('conc1').setIndex('species', true, 0);
            model.component('comp1').physics('tds').feature('conc1').setIndex('c0', 'Ce_0', 0);
            model.component('comp1').physics('tds').feature('conc2').setIndex('species', true, 0);
        end
        
        function model = add_tortuosity_study(model)
            model.study.create('std2');
            model.study('std2').label('Study 2 - Tortuosity');
            
            model.study('std2').create('stat', 'Stationary');
            model.study('std2').feature('stat').set('activate', {'C_e' 'off' ...
                'C_s' 'off' 'phi_l' 'off' 'phi_s' 'off' 'ev' 'off'  ...
                'tds' 'on' 'frame:spatial1' 'on' 'frame:material1' 'on'});
            
            model.sol.create('sol3');
            model.sol('sol3').study('std2');
        end
        
        function model = run_tortuosity_study(model)
            model.study("std2").run();
        end
        
        function [flux, model] = flux_for_tortuosity(model)
            model.result.numerical.create('av1', 'AvLine');
            model.result.numerical('av1').set('intsurface', true);
            model.result.numerical('av1').set('data', 'dset2');
            model.result.numerical('av1').selection.named('geom1_boxsel3');
            model.result.numerical('av1').label('Tortuosity Measurement - Flux');
            model.result.numerical('av1').set('expr', {'tds.tflux_cx'});
            model.result.numerical('av1').set('descr', {'Total flux, x component'});
            model.result.numerical('av1').set('unit', {'mol/(m^2*s)'});
            model.result.table.create('tbl2', 'Table');
            model.result.table('tbl2').comments('Tortuosity Measurement - Flux {av1}');
            model.result.numerical('av1').set('table', 'tbl2');
            model.result.numerical('av1').setResult;
            
            flux = str2double(model.result.table('tbl2').getTableData(false));
            
            % Remove table and line average integral
            model.result.table.remove('tbl2');
            model.result().numerical.remove('av1');

        end
        
        function model = export_electrochem_data(model, num, C, rootpath)
            % export_electrochem_data saves the electrochemical (cell)
            % operation data for microstructure (num) and saves it to a
            % specified filepath as a CSV file.
            %
            % model: COMSOL model
            % num: the microstructure created at "num" iteration
            % C: the C-rate the cell was discharged at
            % rootpath: root of where to save the created CSV file. For
            % instance, "/ROOT/results/%d/electro.csv".
            model.result.export.create('data1', 'Data');
            model.result.export('data1').set('expr', {});
            model.result.export('data1').set('descr', {});
            
            % Specify State-of-Lithiation here
            model.result.export('data1').setIndex('expr', 'C_s/Cs_max * 100', 0);
            model.result.export('data1').setIndex('descr', 'State of Lithiation', 0);
            
            % Make directory if not exist, otherwise, save data
            mkdir(sprintf('%s/results/%d', rootpath, num));
            model.result.export('data1').set('filename', ...
                sprintf('%s/results/%d/electro_%s.csv', ...
                rootpath, num, num2str(C, 2)));
            
            model.result.export('data1').run;
            model.result.export.remove('data1');
        end
    end
end