% This script calibrates the Swisstom measurement data using the assumption
% that the homogenous tank has no frequency differential properties. Thus,
% the measurements taken at the homogenous tank are used as calibration
% data. After calibration an fd reconstruction is performed and the
% resulting fdEIT-image as well as the calibrated voltages are displayed.
%
% All relevant parameters can be specified below.


%% Set parameters

el_no = 32; % Number of measurement electrodes used
nSkip = 0; % Specify the skip pattern
veggie_type = 13; % Specify if only the potato, only the pumpkin or both shall be reconstructed (11, 12 or 13)
inhom_filename = '25.mat'; % Specify the filename of the inhomogenous measurement to be reconstructed
output_name = 'COARSE_adj_78kHz_195kHz_allHPs_homocalibrated_Middle'; % Specify the name of the output file to be saved as

td_reconstruct = true; % Specify if a time differential reconstruction before and after calibration shall be performed
fd_reconstruct = true; % Specify if a frequency differential reconstruction after calibration shall be performed
mesh_size = 'coarse'; % Specify the mesh density ('fine', 'coarse')

% Specify the hyperparameter(s) for the fd-reconstruction
hp_param=[1e-2, 5e-2, 1e-1, 5e-1:1e-1:1, 1.5, 2, 2.5, 3.5, 4, 5, 10, 20, 50];
% hp_param=[1e-2:1e-2:9e-2, 1e-1:1e-1:9e-1, 1:1:10];
% hp_param=[1e-1, 5e-1, 9e-1];
% hp_param=[5e-1];

% Specify the mixed-minimization-parameter(s) for the fd-reconstruction
% beta=[0,0.25,0.5,0.75,1,1.25,2.5,5,7.5,10];
% beta=[0.25,0.5,0.75];
beta=[0.5];

init_iter = 0.0112; % initial conductivity


%% Load data

% Create the measurement path (pattern)
if nSkip == 0
    path = '.\7_Measurement_Evaluation\2020_03_05_swisstom\indexed\adj\';
elseif nSkip == 4
    path = '.\7_Measurement_Evaluation\2020_03_05_swisstom\indexed\sk44\';
else
    error('Please choose a valid measurement pattern.');
end

% Create the measurement path (veggie type)
if veggie_type == 11
    path = strcat(path, '11\');
elseif veggie_type == 12
    path = strcat(path, '13\');
elseif veggie_type == 13
    path = strcat(path, '13\');
else
    error('Please choose a valid vegetable type.');
end

% Create the stimulation pattern
[stim, msel] = mk_stim_patterns(el_no, 1, [0,nSkip+1], [0,nSkip+1], {'no_rotate_meas'}, 1);

% Load the data
load(strcat(path, '20.mat'));
hom_raw = [mean(kHz78, 2), mean(kHz195, 2)];
load(strcat(path, inhom_filename));
inhom_raw = [mean(kHz78, 2), mean(kHz195, 2)];

inhom_raw = [inhom_raw(msel, 1), inhom_raw(msel, 2)]; 
hom_raw = [hom_raw(msel, 1), hom_raw(msel, 2)];

clear kHz78
clear kHz195

%% Calibrate the data
% The data is calibrated under the assumption that the homogenous tank
% should produce exactly the same measurements for both frequencies.

calib = hom_raw(:, 1)./hom_raw(:, 2);
hom_calib = [hom_raw(:, 1), hom_raw(:, 2).*calib];

inhom_calib = [inhom_raw(:, 1), inhom_raw(:, 2).*calib];

%% time-differential reconstruction for testing if everything worked

if td_reconstruct
    
    % Create a fmdl with the dimensions of the tank
    height = 1;
    radius = 0.38;
    el_radius = 0.005;
    el_height = 0.6;
    cyl_shape = [height, radius];
    
    fmdl_recon = ng_mk_cyl_models(cyl_shape, [el_no, el_height], [el_radius, 0, 0]);
    
    fmdl_recon.stimulation = stim;
    
    % Set reconstruction options for GREIT
    opt.imgsz = [64 64];
    opt.Nsim = 500;
    opt.target_size = 0.03;
    opt.target_offset = 0;
    opt.noise_figure = 0.5;
    
    greit_mdl = mk_GREIT_model(fmdl_recon, 0.2, 0.5, opt);
    
    % Reconstruction
    img_solve_raw_f1 = inv_solve(greit_mdl, inhom_raw(:, 1), hom_raw(:, 1));
    img_solve_raw_f2 = inv_solve(greit_mdl, inhom_raw(:, 2), hom_raw(:, 2));
    
    img_solve_f1 = inv_solve(greit_mdl, inhom_calib(:, 1), hom_calib(:, 1));
    img_solve_f2 = inv_solve(greit_mdl, inhom_calib(:, 2), hom_calib(:, 2));
    
    % Take abs of reconstructed data
    img_solve_raw_f1.elem_data = abs(img_solve_raw_f1.elem_data);
    img_solve_raw_f2.elem_data = abs(img_solve_raw_f2.elem_data);
    img_solve_f1.elem_data = abs(img_solve_raw_f1.elem_data);
    img_solve_f2.elem_data = abs(img_solve_raw_f2.elem_data);

    % Display reconstruction
    figure;
    show_slices(img_solve_raw_f1);
    title('td Rekonstruktion 78 kHz');
    
    figure;
    show_slices(img_solve_raw_f2);
    title('td Rekonstruktion 195 kHz');
    
    figure;
    show_slices(img_solve_f1);
    title('Kalibrierte td Rekonstruktion 78 kHz');
    
    figure;
    show_slices(img_solve_f2);
    title('Kalibrierte td Rekonstruktion 195 kHz');    
end

%% frequency differential reconstruction using calibrated data

% Create fmdl for fd-reconstruction
height = 0;
radius = 14.5;
el_radius = 0.5;
el_height = 0;

% Define mesh density
if strcmp(mesh_size, 'coarse')
    max_mesh_size = 2;
    max_mesh_size_electrodes = 0.5;
elseif strcmp(mesh_size, 'fine')
    max_mesh_size = 1;
    max_mesh_size_electrodes = 0.05;
end

cyl_shape = [height, radius, max_mesh_size];
fmdl_recon = ng_mk_cyl_models(cyl_shape, [el_no, el_height], [el_radius, 0, max_mesh_size_electrodes]);

fmdl_recon.stimulation = stim;

%% Reconstruct

% Reconstruct for all lambdas and betas specified above
if fd_reconstruct
        for jj = 1:length(beta)
            for ii = 1:length(hp_param)
                
                imdl = my_weighted_invprob_properties(fmdl_recon, init_iter, hp_param(ii), true,beta(jj));
                img_solve(ii)= inv_solve(imdl, inhom_calib(:, 1), inhom_calib(:, 2));
                
                disp('o-----------------------------------------------------------------------------------o')
                disp('Finished reconstructin for hyperparameter: ')
                disp(num2str(ii));
                disp('Last reconstruction used beta number: ')
                disp(num2str(jj));
                disp('o-----------------------------------------------------------------------------------o')
            end
            
            % Save the reconstrucion
            save(strcat(output_name,'.mat'),'img_solve','beta','hp_param');
            
        end
        
        figure;
        show_slices(img_solve);        
end


%% Display the voltages before and after calibration

figure;

subplot(2,2,1)
plot(abs(hom_raw(:, 2)))
ylim([0 3.5e-4]);
title('Homogen, unkalibriert, 195kHz');

subplot(2,2,2)
plot(abs(hom_calib(:, 2)))
ylim([0 3.5e-4]);
title('Homogen, kalibriert, 195kHz');

subplot(2,2,3)
plot(abs(inhom_raw(:, 2)))
ylim([0 3.5e-4]);
title('Inhomogen, unkalibriert, 195kHz');

subplot(2,2,4)
plot(abs(inhom_calib(:, 2)))
ylim([0 3.5e-4]);
title('Inhomogen, kalibriert, 195kHz');

figure;
subplot(2,1,1)
plot(abs(inhom_raw(:, 2)-inhom_raw(:, 1)))
title('Inhomogen, Differenz, unkalibriert');

subplot(2,1,2)
plot(abs(inhom_calib(:, 2)-inhom_calib(:, 1)))
title('Inhomogen, Differenz, kalibriert');