% This script will reconstruct measurements taken with the AixTOM
% measurement system in a time differential way. All parameters for the
% reconstruction as well as the measurement data can be specified by the
% user at the beginning of the script.

%% Clear relevant variables in the workspace
clear voltages;
clear filenames;
clear measurement;
clear measurement_data;
clear measurement_title;
close all

%% Set the reconstruction parameters

el_no = 16; % Number of electrodes used for measuring
nSkip = 4; % Specify the skip pattern used
frequency = 2; % Specify which of the 2 frequecies that were used for measuring shall be used for td reconstruction (1 or 2)
output_name = 'td_recon_sk44_newtank_13_72016_'; % Name of the file that will be stored
path = '.\7_Measurement_Evaluation\Older\2020_02_20\Indexed\Sk44\New_Tank\12\'; % Measurement path

take_mean = true; % Shall a mean of all frames in this measurement be taken?

% Specify if a conintuous set of measurement data shall be loaded (e.g. all
% files in a folder like 20, 21, 22, 23, 24, 25). This required that the
% indexing of the data was done correctily since the first measurement will
% always automatically be interpreted as the homogenous measurement.
continuous = false; 

% Mark Start and end filename here, if continuous is true. i_start must
% always be the homogenous measurement
i_start = 20;
i_finish = 25;

% Specify filenames for the homogenous and the inhomogenous measurements
% here if continuous is false
homogenous_filename = '20';
filenames = {'21','22','23','24'};


%% Load the data

if continuous
    for i = i_start:i_finish
        
        % The measurements are stored in the measruement cell-array. The
        % indexing is chosen such that the first element is always the
        % homogenous file. For this to work the homogenous file must be the
        % first measurement specified in for the continuous loading above.
        measurement{i+1-i_start} = AixTOM_readFile(strcat(path, num2str(i), '.txt'));
        measurement_data{i+1-i_start} = load(strcat(path, num2str(i), '.mat'));

        meas_title{i+1-i_start} = measurement_data{i+1-i_start}.deviceData.comment;
        
        % Take mean if desired.
        if take_mean
            voltages{i+1-i_start} = mean(measurement{i+1-i_start}(:, :, frequency), 1, 'omitnan')';
        else
            voltages{i+1-i_start} = measurement{i+1-i_start}(:, :, frequency)';
        end
        
    end
else
    
    % First, the homogenous measurement is loaded and stored as the first
    % item in the cell array
    measurement{1} = AixTOM_readFile(strcat(path, homogenous_filename, '.txt'));
    measurement_data{1} = load(strcat(path, homogenous_filename, '.mat'));
    meas_title{1} = measurement_data{1}.deviceData.comment;
    
    if take_mean
        voltages{1} = mean(measurement{1}(:, :, frequency), 1, 'omitnan')';
    else
        voltages{1} = measurement{1}(:, :, frequency)';
    end
    
    % Now all the other files are loaded
    for i = 2:(length(filenames) + 1)
        
        measurement{i} = AixTOM_readFile(strcat(path, filenames{i-1}, '.txt'));
        measurement_data{i} = load(strcat(path, filenames{i-1}, '.mat'));
        meas_title{i} = measurement_data{i}.deviceData.comment;
        
        if take_mean
            voltages{i} = mean(measurement{i}(:, :, frequency), 1, 'omitnan')';
        else
            voltages{i} = measurement{i}(:, :, frequency)';
        end
        
    end
    
end

%% Generate forward model

% Generate a 3D fmdl that has the same proportions as the tank if the tank
% is filled with 25l of water
height = 1;
radius = 0.38;
el_radius = 0.005;
el_height = 0.6;
cyl_shape = [height, radius];

fmdl_recon = ng_mk_cyl_models(cyl_shape, [el_no, el_height], [el_radius, 0, 0]);

% Generate adjacent stimulation pattern
fmdl_recon.stimulation = mk_stim_patterns(el_no ,1,[15 nSkip],[15 nSkip],{'no_rotate_meas'},1);

%% Reconstruct all measurements

for k = 2:length(measurement)
    
    % Set reconstruction options for GREIT
    opt.imgsz = [64 64];
    opt.Nsim = 500;
    opt.target_size = 0.03;
    opt.target_offset = 0;
    opt.noise_figure = 0.5;
    
    % Make a GREIT model for reconstruction
    greit_mdl = mk_GREIT_model(fmdl_recon, 0.2, 0.5, opt);
    
    % Reconstruction
    img_solve = inv_solve(greit_mdl, voltages{k}, voltages{1});
    current_title = meas_title{k};
    
    disp('o===================================================================================o')
    disp('Finished reconstruction of measurement: ')
    disp(num2str(k + 19))
    disp('o===================================================================================o')
    
    % Display the reconstruction result
    figure;
    img_solve.elem_data = abs(img_solve.elem_data);
    show_slices(img_solve);
    eidors_colourbar(img_solve);
    title(meas_title{k}, 'Interpreter', 'none');
    
    % Display the differential measurement voltages
    figure;
    plot(abs(voltages{k}));
    title(strcat(meas_title{k}, ', Spannungen (abs)'), 'Interpreter', 'none');

end

