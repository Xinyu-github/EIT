% This script will reconstruct measurements taken with the SwissTOM EIT
% Pioneer Set in a time differential way. All parameters for the
% reconstruction as well as the measurement data can be specified by the
% user at the beginning of the script.

%% Set Parameters

% Set measurement Path
path = 'C:\Code\fdEIT\7_Measurement_Evaluation\2020_03_05_swisstom\indexed\sk44\11\';

el_no = 32; % Number of measurement electroded
nSkip = 4; % Specify the skip pattern

%% Load the data

load(strcat(path, '20.mat'));
hom_low = mean(kHz78, 2);
hom_high = mean(kHz195, 2);

load(strcat(path, '21.mat'));
inhom_low = mean(kHz78, 2);
inhom_high = mean(kHz195, 2);

%% Reconstruction

% Create a fmdl with the tank dimensions
height = 1;
radius = 0.38;
el_radius = 0.005;
el_height = 0.6;
cyl_shape = [height, radius];

fmdl_recon = ng_mk_cyl_models(cyl_shape, [el_no, el_height], [el_radius, 0, 0]);

[stim, msel] = mk_stim_patterns(el_no, 1, [0,nSkip+1], [0,nSkip+1], {'no_rotate_meas'}, 1);
fmdl_recon.stimulation = stim;

% Set reconstruction options for GREIT
opt.imgsz = [64 64];
opt.Nsim = 500;
opt.target_size = 0.03;
opt.target_offset = 0;
opt.noise_figure = 0.5;

greit_mdl = mk_GREIT_model(fmdl_recon, 0.2, 0.5, opt);

% Reconstruction
img_solve = inv_solve(greit_mdl, inhom_high(msel), hom_high(msel));

% Display the result
figure;
img_solve.elem_data = abs(img_solve.elem_data);
show_slices(img_solve);
img_solve.calc_colours.ref_level =  [];
img_solve.calc_colours.clim      =  [];

% Display the voltages
figure;
plot(abs(inhom_high(msel)));
title('Kartoffel, 78kHz, El_1, Spannungen (abs)', 'Interpreter', 'none');
saveas(gcf, 'voltages_11_21_195kHz_sk44_swisstom.svg');
xlim([1 928])



