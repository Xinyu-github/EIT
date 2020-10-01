% This script can be used to reconstruct data taken with the Dräger EEK2
% measurement device. The reconstruction is done time differential using
% the GREIT algorithm. 
%
% Please notice that the EEK2 data is taken in a time continuous way. Thus,
% the different vegetable positions have to be separated from each other
% manually. This has already been done in this script for the existing
% measurement. However the separation has to be done again if new data is
% measured.

clear
close all

%% Set parameters

% Measurement path and filename
path = '.\7_Measurement_Evaluation\Older\2020_03_02-03\Not Indexed\EEK_2\Adj\New_Tank\11\';
filename = '11_eik2_20200302_adj_potato_01_001_me01_eiu.eit';

el_no = 16; % Specify electrode number


%% Load data

[data, auxdata, stim] = eidors_readdata(strcat(path, filename), 'DRAEGER-EIT', []);

% Separate the different measurements in the time-continuous measurement
homogen = mean(data(:, 1:1848), 2);
mitte = mean(data(:, 2314:4620), 2);
el_1 = mean(data(:, 8040:9091), 2);
el_5 = mean(data(:, 12360:13660), 2);
el_9 = mean(data(:, 16820:18460), 2);
el_13 = mean(data(:, 22310:end), 2);


% Plot the different measured voltages

% figure;
% subplot(2, 3, 1);
% plot(homogen);
% title('Homogen');
% xlabel('Measurement Number');
% ylabel('Measured Voltage [V]');
% 
% subplot(2, 3, 2);
% plot(mitte);
% title('Mitte');
% xlabel('Measurement Number');
% ylabel('Measured Voltage [V]');
% 
% subplot(2, 3, 3);
% plot(el_1);
% title('Elektrode 1');
% xlabel('Measurement Number');
% ylabel('Measured Voltage [V]');
% 
% subplot(2, 3, 4);
% plot(el_5);
% title('Elektrode 5');
% xlabel('Measurement Number');
% ylabel('Measured Voltage [V]');
% 
% subplot(2, 3, 5);
% plot(el_9);
% title('Elektrode 9');
% xlabel('Measurement Number');
% ylabel('Measured Voltage [V]');
% 
% subplot(2, 3, 6);
% plot(el_13);
% title('Elektrode 13');
% xlabel('Measurement Number');
% ylabel('Measured Voltage [V]');

%% Create forward model

% Generate a 3D fmdl that has the same proportions as the tank if the tank
% is filled with 25l of water
height = 1;
radius = 0.38;
el_radius = 0.005;
el_height = 0.6;
cyl_shape = [height, radius];

% Create stimulation pattern
fmdl_recon = ng_mk_cyl_models(cyl_shape, [el_no, el_height], [el_radius, 0, 0]);
fmdl_recon.stimulation = stim;    


%% Reconstruct

% Set reconstruction options for GREIT
opt.imgsz = [64 64];
opt.Nsim = 500;
opt.target_size = 0.03;
opt.target_offset = 0;
opt.noise_figure = 0.5;

greit_mdl = mk_GREIT_model(fmdl_recon, 0.2, 0.5, opt);

% Reconstruction
img_solve_mitte = inv_solve(greit_mdl, mitte, homogen);
img_solve_el1 = inv_solve(greit_mdl, el_1, homogen);
img_solve_el5 = inv_solve(greit_mdl, el_5, homogen);
img_solve_el9 = inv_solve(greit_mdl, el_9, homogen);
img_solve_el13 = inv_solve(greit_mdl, el_13, homogen);

% Display the reconstructions
figure;
show_slices(img_solve_el1);
title('el 1');
figure;
show_slices(img_solve_el5);
title('el 5');
figure;
show_slices(img_solve_el9);
title('el 9');
figure;
show_slices(img_solve_el13);
title('el 13');
figure;
show_slices(img_solve_mitte);
title('mitte');


