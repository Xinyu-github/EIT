% This script will reconstruct measurements taken with the SwissTOM EIT
% Pioneer Set in a frequency differential way. All parameters for the
% reconstruction as well as the measurement data can be specified by the
% user at the beginning of the script.

%% Clear relevant variables in the workspace
clear;
close all;

%% Set reconstruction parameters

display_fem = false; % Specify if the FEM used as a forward model shall be displayed
reconstruct = true; % Specify if an fd reconstruction shall be performed
td_reconstruct = false; % Specify if a td reconstruction of the measurement data shall be performed
calibrate = true; % Specify if the raw measurement data shall be calculated using the measurements of the calibration network developed at MEDIT
display_calibration_voltages = false; % specify if the voltages before and after shall be displayed

% Specify if a calibration using the homogenous tank as a calibrator shall 
% be performed on top of the calibration done using the clibration network.
% This parameter can only be true if the parameter 'calibrate' is also
% true.
hom_calibrate = false; 

el_no = 32; % Number of measurement electrodes
nSkip = 0; % Specify the skip pattern
mesh_size = 'coarse'; % Define the density of the mesh ('fine' or 'coarse')
normalize = false; % Specify if the data shall be normalized

% Set the name for the export-file
output_name = 'COARSE_sk44_78kHz_195kHz_allHPs_';

% Specify path and filename of the calibration measurements
cali_path = 'C:\Code\fdEIT\7_Measurement_Evaluation\2020_03_05_swisstom\indexed\adj\calibration\';
cali_filename_8Ohm = 'calibrator_high_ohm.mat';
cali_filename_0_1Ohm = 'calibrator_low_ohm.mat';

% Specify the measurement directory and the names of the files
path = 'C:\Code\fdEIT\7_Measurement_Evaluation\2020_03_05_swisstom\indexed\adj\13\';
filenames = {'25'};

% Set the reconstruction parameters
% Select variant hyper-parameters and different absolute-mixed parameters
% hp_param_all = exp((1:2:40))/(10^7); hp_param = hp_param_all(1:12);
% hp_param=[1e-2, 5e-2, 1e-1, 5e-1:1e-1:1, 1.5, 2, 2.5, 3.5, 4, 5, 10, 20, 50];
% hp_param=[1e-2:1e-2:9e-2, 1e-1:1e-1:9e-1, 1:1:10];
% hp_param=[1e-1, 5e-1, 9e-1];
hp_param=[5e-1];
% beta=[0,0.25,0.5,0.75,1,1.25,2.5,5,7.5,10];
% beta=[0.25,0.5,0.75];
beta=[0.5];

init_iter = 0.0112; % initial conductivity

%% Load data and generate stimulation pattern

% Generate stimulation pattern
[stim, msel] = mk_stim_patterns(el_no, 1, [0,nSkip+1], [0,nSkip+1], {'no_rotate_meas'}, 1);

% Load measurement data and set the names of the individual measurements.
% If you wish to add data that doesn't follow the indexing convention, the
% names will be wrong!
for i = 1:length(filenames)
    
    measurement{i} = load(strcat(path, filenames{i}, '.mat'));
    
    switch filenames{i}
        case '20'
            voltages1{i} = mean(measurement{i}.kHz78, 2, 'omitnan')';
            voltages2{i} = mean(measurement{i}.kHz195, 2, 'omitnan')';
            titles{i} = 'Homogenous';
        case '21'
            voltages1{i} = mean(measurement{i}.kHz78, 2, 'omitnan')';
            voltages2{i} = mean(measurement{i}.kHz195, 2, 'omitnan')';
            titles{i} = 'Electrode_1';
        case '22'
            voltages1{i} = mean(measurement{i}.kHz78, 2, 'omitnan')';
            voltages2{i} = mean(measurement{i}.kHz195, 2, 'omitnan')';
            titles{i} = 'Electrode_9';
        case '23'
            voltages1{i} = mean(measurement{i}.kHz78, 2, 'omitnan')';
            voltages2{i} = mean(measurement{i}.kHz195, 2, 'omitnan')';
            titles{i} = 'Electrode_17';
        case '24'
            voltages1{i} = mean(measurement{i}.kHz78, 2, 'omitnan')';
            voltages2{i} = mean(measurement{i}.kHz195, 2, 'omitnan')';
            titles{i} = 'Electrode_25';
        case '25'
            voltages1{i} = mean(measurement{i}.kHz78, 2, 'omitnan')';
            voltages2{i} = mean(measurement{i}.kHz195, 2, 'omitnan')';
            titles{i} = 'Middle';
    end        
 
end

%% Calibration

if calibrate
    
    % load calibration data
    cali_meas_8Ohm_raw = load(strcat(cali_path, cali_filename_8Ohm));
    cali_meas_0_1Ohm_raw = load(strcat(cali_path, cali_filename_0_1Ohm));
    
    cali_meas_8Ohm = [mean(cali_meas_8Ohm_raw.kHz78, 2), mean(cali_meas_8Ohm_raw.kHz195, 2)];
    cali_meas_0_1Ohm = [mean(cali_meas_0_1Ohm_raw.kHz78, 2), mean(cali_meas_0_1Ohm_raw.kHz195, 2)];
    cali_meas_8Ohm = [cali_meas_8Ohm(msel, 1), cali_meas_8Ohm(msel, 2)];
    cali_meas_0_1Ohm = [cali_meas_0_1Ohm(msel, 1), cali_meas_0_1Ohm(msel, 2)];
    
    % Calibrate the data using the designated calibration function
    [calib_f1, calib_f2, calib_f1_pol2cart, calib_f2_pol2cart] = swissTOM_calibrate(cali_meas_8Ohm , cali_meas_0_1Ohm, (voltages1{1}(msel))', (voltages2{1}(msel))', display_calibration_voltages, 'pumpkin and potato');
    
    % Calibrate the homogenous measurement for the td-reconstruction
    homogen = load(strcat(path, '20.mat'));
    mean_homogen_f1 = mean(homogen.kHz78, 2);
    mean_homogen_f2 = mean(homogen.kHz195, 2);
    mean_homogen_f1 = mean_homogen_f1(msel);
    mean_homogen_f2 = mean_homogen_f2(msel);
    [calib_hom_f1, calib_hom_f2, calib_hom_f1_pol2cart, calib_hom_f2_pol2cart] = swissTOM_calibrate(cali_meas_8Ohm , cali_meas_0_1Ohm, mean_homogen_f1, mean_homogen_f2, display_calibration_voltages, 'homogen');

    if hom_calibrate
        
        figure;
        subplot(2,1,1)
        plot(abs(calib_f2_pol2cart-calib_f1_pol2cart));
        title('abs of difference before homocal');
        
        % Perform a homogenous calibration on top of the calibration that
        % was already performed
        hom_calib = calib_hom_f1_pol2cart./calib_hom_f2_pol2cart;
        calib_hom_f2_pol2cart = calib_hom_f2_pol2cart.*hom_calib;
        calib_f2_pol2cart = calib_f2_pol2cart.*hom_calib;

        subplot(2,1,2)
        plot(abs(calib_f2_pol2cart-calib_f1_pol2cart));
        title('abs of difference after homocal');
        
    end

end


%% Reconstruct td to check if calibration was successfull

if td_reconstruct
    
    % Create fmdl with the dimenstions of the measurement tank
    height = 1;
    radius = 0.38;
    el_radius = 0.005;
    el_height = 0.6;
    cyl_shape = [height, radius];
    
    fmdl_recon = ng_mk_cyl_models(cyl_shape, [32, el_height], [el_radius, 0, 0]);
    
    % Normalise if desired
    if normalize
        fmdl_recon = mdl_normalize(fmdl_recon,  1);
    end
    
    fmdl_recon.stimulation = stim;
    
    % Set reconstruction options for GREIT
    opt.imgsz = [64 64];
    opt.Nsim = 500;
    opt.target_size = 0.03;
    opt.target_offset = 0;
    opt.noise_figure = 0.5;
    
    greit_mdl = mk_GREIT_model(fmdl_recon, 0.2, 0.5, opt);
    
    % Reconstruction
    img_solve_mitte_cal = inv_solve(greit_mdl, calib_f1_pol2cart, calib_hom_f1_pol2cart);
    
    figure;
    img_solve_mitte_cal.elem_data = abs(img_solve_mitte_cal.elem_data);
    show_slices(img_solve_mitte_cal);
    eidors_colourbar(img_solve_mitte_cal);
    title('mitte, cal');
    
    % Reconstruction
    img_solve_mitte = inv_solve(greit_mdl, (voltages1{1}(msel))', mean_homogen_f1);
    
    figure;
    img_solve_mitte.elem_data = abs(img_solve_mitte.elem_data);
    show_slices(img_solve_mitte);
    eidors_colourbar(img_solve_mitte);
    title('mitte, no cal');
end

%% Generate forward model for fd reconstruction

% Create fmdl according to the tank dimensions
height = 0;
radius = 14.5;
el_radius = 0.5;
el_height = 0;

% Set the mesh density
if strcmp(mesh_size, 'coarse')
    max_mesh_size = 2;
    max_mesh_size_electrodes = 0.5;
elseif strcmp(mesh_size, 'fine')
    max_mesh_size = 1;
    max_mesh_size_electrodes = 0.05;
end

cyl_shape = [height, radius, max_mesh_size];
fmdl_recon = ng_mk_cyl_models(cyl_shape, [el_no, el_height], [el_radius, 0, max_mesh_size_electrodes]);

% Display FEM if desired
if display_fem
    figure;
    show_fem(fmdl_recon, [0 1 0]);
    axis off
end

if normalize
    fmdl_recon = mdl_normalize(fmdl_recon,  1);
end

fmdl_recon.stimulation = stim;

%% Frequency differential reconstruction

if reconstruct
    
    % Reconsruct all measurements
    for k = 1:length(measurement)
        
        % Define whether to use the calibrated or non-calibrated voltages
        if calibrate
            high_volt = calib_f2_pol2cart;
            low_volt = calib_f1_pol2cart;
        else
            high_volt = (voltages2{1}(msel))';
            low_volt = (voltages1{1}(msel))';
        end
        
        % Reconstruct for all lambdas and betas
        for jj = 1:length(beta)
            for ii = 1:length(hp_param)
                
                imdl = my_weighted_invprob_properties(fmdl_recon, init_iter, hp_param(ii), true,beta(jj));
                img_solve(ii)= inv_solve(imdl, low_volt, high_volt);
                
                disp('o-----------------------------------------------------------------------------------o')
                disp('Finished reconstructin for hyperparameter: ')
                disp(num2str(ii));
                disp('Last reconstruction used beta number: ')
                disp(num2str(jj));
                disp('o-----------------------------------------------------------------------------------o')
            end
            
            current_title = titles{k};
            save(strcat(output_name, titles{k},'.mat'),'img_solve','beta','hp_param', 'current_title');
            
        end
        
        disp('o===================================================================================o')
        disp('Finished reconstruction of measurement: ')
        disp(k)
        disp('o===================================================================================o')
        
        figure;
        img_solve.elem_data = abs(img_solve.elem_data);
        show_slices(img_solve);
        eidors_colourbar(img_solve);
        title(titles{k});
        
        beep
    end

end
