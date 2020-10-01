% This script will reconstruct measurements taken with the AixTOM
% measurement system in a frequency differential way. All parameters for the
% reconstruction as well as the measurement data can be specified by the
% user at the beginning of the script.

%% Clear relevant variables in the workspace
clear voltages1;
clear voltages2;
clear filenames;
clear measurement;
clear measurement_data;
clear measurement_title;

%% Set the reconstruction parameters

el_no = 16; % Number of electrodes used for the measurement
nSkip = 4; % Specify the Skip-Pattern that was used
mesh_size = 'coarse'; % Define the density of the mesh ('fine' or 'coarse')
output_name = 'fd_recon_sk44_newtank_11_21_TEST'; % Specify the name for the export file
path = './7_Measurement_Evaluation/Older/2020_02_20/Indexed/Sk44/New_Tank/12/'; % Specify the measurement path

take_mean = true; % Specify if a mean of all measurement frames shall be taken
normalize = false; % Specify if the forward model shall be normalised

% Specify if a conintuous set of measurement data shall be loaded (e.g. all
% files in a folder like 21, 22, 23, 24, 25).
continuous = false;

% Mark Start and end filename here, if continuous is true
i_start = 21;
i_finish = 25;

% Specify filenames of the measurements to be reconstructed if continuous 
% is false
filenames = {'20', '21'};

% Specify the reconstruction parameters
% hp_param_all = exp((1:2:40))/(10^7); hp_param = hp_param_all(1:12);
% hp_param=[1e-2, 5e-2, 1e-1, 5e-1:1e-1:1, 1.5, 2, 2.5, 3.5, 4, 5, 10, 20, 50];
% hp_param=[1e-2:1e-2:9e-2, 1e-1:1e-1:9e-1, 1:1:10];
% hp_param=[1e-1, 5e-1:1e-1:1, 20];
hp_param=[5e-1];
% beta=[0,0.25,0.5,0.75,1,1.25,2.5,5,7.5,10];
% beta=[0.25,0.5,0.75];
beta=[0.5];

init_iter = 0.0112; % initial conductivity


%% Load the specified measurement data

if continuous
    for i = i_start:i_finish
        
        % The measurements are stored in the measurement cell-array such
        % that the first of the continous filenames is the first item in
        % the array
        measurement{i+1-i_start} = AixTOM_readFile(strcat(path, num2str(i), '.txt'));
        measurement_data{i+1-i_start} = load(strcat(path, num2str(i), '.mat'));

        meas_title{i+1-i_start} = measurement_data{i+1-i_start}.deviceData.comment;
        
        % Take mean if desired
        if take_mean
            voltages1{i+1-i_start} = mean(measurement{i+1-i_start}(:, :, 1), 1, 'omitnan')';
            voltages2{i+1-i_start} = mean(measurement{i+1-i_start}(:, :, 2), 1, 'omitnan')';
        else
            voltages1{i+1-i_start} = measurement{i+1-i_start}(:, :, 1)';
            voltages2{i+1-i_start} = measurement{i+1-i_start}(:, :, 2)';
        end
        
    end
else
    
    for i = 1:length(filenames)
        
        % If continuous if false, load only the specified files
        measurement{i} = AixTOM_readFile(strcat(path, filenames{i}, '.txt'));
        measurement_data{i} = load(strcat(path, filenames{i}, '.mat'));
        meas_title{i} = measurement_data{i}.deviceData.comment;
        
        if take_mean
            voltages1{i} = mean(measurement{i}(:, :, 1), 1, 'omitnan')';
            voltages2{i} = mean(measurement{i}(:, :, 2), 1, 'omitnan')';
        else
            voltages1{i} = measurement{i}(:, :, 1)';
            voltages2{i} = measurement{i}(:, :, 2)';
        end
        
    end
    
end


%% Generate forward model

% Generate a forward model with the radius of the tank in cm
height = 0;
radius = 14.5;
el_radius = 0;
el_height = 0;
max_mesh_size = 1;
cyl_shape = [height, radius, max_mesh_size];

% Define the mesh density
if strcmp(mesh_size, 'coarse')
    max_mesh_size = 2;
    max_mesh_size_electrodes = 0.5;
elseif strcmp(mesh_size, 'fine')
    max_mesh_size = 1;
    max_mesh_size_electrodes = 0.05;
end

cyl_shape = [height, radius, max_mesh_size];
fmdl_recon = ng_mk_cyl_models(cyl_shape, [el_no, el_height], [el_radius, 0, max_mesh_size_electrodes]);

figure;
show_fem(fmdl_recon);

% normalise if desired
if normalize
    fmdl_recon = mdl_normalize(fmdl_recon,  1);
end

% Generate adjacent stimulation pattern
fmdl_recon.stimulation = mk_stim_patterns(el_no ,1,[15 nSkip],[15 nSkip],{'no_rotate_meas'},1);



%% Reconstruct

% Reconstruct all measurements stored in the measurement cell array
for k = 1:length(measurement)
    
    % Reconstruct for all betas and lambdas specified above
    for jj = 1:length(beta)
        for ii = 1:length(hp_param)
            
            imdl    = my_weighted_invprob_properties(fmdl_recon, init_iter, hp_param(ii),...
                true,beta(jj));
            img_solve(ii)= inv_solve(imdl, voltages1{k}, voltages2{k});
            
            % Display which hyperparameter and bete were just reconstructed
            disp('o-----------------------------------------------------------------------------------o')
            disp('Finished reconstructin for hyperparameter: ')
            disp(num2str(ii))
            disp('Last reconstruction used beta number: ')
            disp(num2str(jj))
            disp('o-----------------------------------------------------------------------------------o')
        end
        
        % Save the reconstruction
        if continuous
            current_title = meas_title{k};
            save(strcat(output_name ,num2str(k + 20),'.mat'),'img_solve','beta','hp_param', 'current_title');
        else
            current_title = meas_title{k};
            save(strcat(output_name, filenames{k},'.mat'),'img_solve','beta','hp_param', 'current_title');
        end
            
    end
    
    disp('o===================================================================================o')
    disp('Finished reconstruction of measurement: ')
    disp(num2str(k + 20))
    disp('o===================================================================================o')
   
    figure;
    show_slices(img_solve);
    title(meas_title{k});
    
    beep
end

