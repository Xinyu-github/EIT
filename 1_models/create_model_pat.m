function [fmdl, img_sim1, img_sim2, img_diff] = create_model_pat(freq_1, freq_2, pathology, elec_num, grid_size, nSkip,ratio)
% This function creates models for testing algorithms on different
% pathologies. The pathologies that can be modelled are specified by the
% input argument 'pathology'. This function is NOT part of the EIDORS
% package and was developed at the MEDIT institute.
%
% Syntax: [fmdl, img_sim1, img_sim2, img_diff] = create_model_pat(freq_1, freq_2, pathology, elec_num, grid_size, nSkip)
%
% Input Parameters:
% freq_1:       Lower frequency to simulate the model pathology at.
% freq_2:       Higher frequency to simulate the mondel pathology at.
% pathology:    This parameter specifies the pathology for which a model is
%               going to be created. The function supports a healthy thorax
%               as well as 9 different pathologies. The parameter can be
%               specified as one of the strings listed below.
%               'healthy':        Healthy Lung
%               'at_left':        Full atelectasis of left lung
%               'at_right':       Full atelectasis of right lung
%               'ed_left':        Edema in left lung
%               'ed_right':       Edema in right lung
%               'ed_low_right':   Edema in lower third of right lung
%               'ed_low_left':    Edema in lower third of left lung
%               'peri_eff':       Pericardial effusion
%               'pleu_eff_left':  Pleural effusion on the left side of the thorax
%               'pleu_eff_right': Pleural effusion on the right side of the thorax
%               'pneu_left':      Pneumothorax on the left side
% elec_num:     Number of electrodes to be placed around the thorax.
% grid_size:    Maximum size of the grid of the fem that is used to model
%               the thorax.
% nSkip:        Number of electrodes between the two stimulation electrodes
%               and each pair of measurement electrodes. This parameter
%               specifies the skip pattern that is used for generating the
%               measurement voltages.
% ratio:        Ratio between lung tissue and extracellular fluid in edema
%               case
% Output Parameters:
% fmdl:         The forward fem that is used to model the specified
%               pathology.
% img_sim1:     EIDORS image with conductivities at lower simulation
%               frequency.
% img_sim2:     EIDORS image with conductivities at higher simulation
%               frequency.
% img_diff:     EIDORS image of the frequency differential conductivities.


%% Load tissue conductivities

tissues_1 = getTissues_mod(freq_1);
tissues_2 = getTissues_mod(freq_2); 

% Get muscle conductivity as background at both frequencies
s_bkgnd_1 = round(tissues_1.Muscle_combined,4,'decimal'); 
s_bkgnd_2 = round(tissues_2.Muscle_combined,4,'decimal'); 

% Get conductivites at frequency 1
s_lung_1 = round(tissues_1.LungInflated_combined,4,'decimal'); 
s_heart_1 = round(tissues_1.Heart_combined,4,'decimal'); 
s_spine_1 = round(tissues_1.BoneAussen_combined,4,'decimal');
s_fluid_1 = round(tissues_1.fluid_combined,4,'decimal'); 
s_blood_1 = round(tissues_1.Blood_combined, 4, 'decimal');

% Get conductivitied at frequency 2
s_lung_2 = round(tissues_2.LungInflated_combined,4,'decimal'); 
s_heart_2 = round(tissues_2.Heart_combined,4,'decimal'); 
s_spine_2 = round(tissues_2.BoneAussen_combined,4,'decimal');
s_fluid_2 = round(tissues_2.fluid_combined,4,'decimal');
s_blood_2 = round(tissues_2.Blood_combined, 4, 'decimal');

try 
    alpha = ratio;
catch
    alpha = 0;
end
%% Make model for the specified pathology

if strcmp(pathology, 'at_left')
    fmdl = create_thorax_fem_simon(elec_num, nSkip, grid_size, 3, 0, 0);
elseif strcmp(pathology, 'at_right')
    fmdl = create_thorax_fem_simon(elec_num, nSkip, grid_size, 2, 0, 0);  
elseif strcmp(pathology, 'pneu_left')
    fmdl = create_thorax_fem_simon(elec_num, nSkip, grid_size, 4, 0, 0); 
else
    fmdl = create_thorax_fem_simon(elec_num, nSkip, grid_size, 1, 0, 0);
end

% Create one image for each frequency and a placeholder for the
% differential image
img_sim1 = mk_image( fmdl, s_bkgnd_1 );  
img_sim2 = mk_image( fmdl, s_bkgnd_2 ); 
img_diff = img_sim1;

%% Create organ shapes

heart_center = [0.12,0.25]; % Assume circular heart; placement by center 
spine_center = [0.01,-0.47]; % Assume ellpictical spine, placement by center
lung_center_right = [0.15,0.1]; % Assume cropped circle as lung, placement by center (left lung is mirrored)

heart_radius = 0.27; % Radius of the circular heart
organ_fn = create_simple_organ_shape( heart_center,heart_radius, spine_center, lung_center_right ); % Load the prespecified organ shapes

% Set the functions that define the organ shapes
fn_heart = organ_fn.heart;
fn_spine = organ_fn.spine;
fn_inner_heart = organ_fn.inner_heart;

%% Apply heart and spine conductivities depending on the pathology
if strcmp(pathology, 'peri_eff')
    img_sim1.elem_data(:) = s_bkgnd_1 ...
        +( (s_blood_1)- (s_bkgnd_1) )*elem_select(img_sim1.fwd_model, fn_heart)...
        +( (s_spine_1)- (s_bkgnd_1) )*elem_select(img_sim1.fwd_model, fn_spine);
    img_sim2.elem_data(:) = s_bkgnd_2 ... 
        +( (s_blood_2)- (s_bkgnd_2) )*elem_select(img_sim2.fwd_model, fn_heart)...
        +( (s_spine_2)- (s_bkgnd_2) )*elem_select(img_sim2.fwd_model, fn_spine);
    img_sim1.elem_data(:) = img_sim1.elem_data(:) + (elem_select(img_sim1.fwd_model, fn_inner_heart)*(s_heart_1 - s_blood_1));
    img_sim2.elem_data(:) = img_sim2.elem_data(:) + (elem_select(img_sim2.fwd_model, fn_inner_heart)*(s_heart_2 - s_blood_2));
    spine_indx = find(round(img_sim1.elem_data(:),4,'decimal')==round(s_spine_1,4,'decimal'));    
else
    img_sim1.elem_data(:) = s_bkgnd_1 ...
        +( (s_heart_1)- (s_bkgnd_1) )*elem_select(img_sim1.fwd_model, fn_heart)...
        +( (s_spine_1)- (s_bkgnd_1) )*elem_select(img_sim1.fwd_model, fn_spine);
    img_sim2.elem_data(:) = s_bkgnd_2 ... 
        +( (s_heart_2)- (s_bkgnd_2) )*elem_select(img_sim2.fwd_model, fn_heart)...
        +( (s_spine_2)- (s_bkgnd_2) )*elem_select(img_sim2.fwd_model, fn_spine);
    spine_indx = find(round(img_sim1.elem_data(:),4,'decimal')==round(s_spine_1,4,'decimal'));
end

%% Apply lung conductivities depending on the pathology

if strcmp(pathology, 'at_left')
    
    % Apply right lung only
    lung_elem_indx = fmdl.mat_idx{1,2};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
elseif strcmp(pathology, 'at_right')
    
    % Apply left lung only
    lung_elem_indx = fmdl.mat_idx{1,2};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
elseif strcmp(pathology, 'ed_left')
    s_mix_1 = conductivity_mixture(s_lung_1,s_fluid_1,alpha);
    s_mix_2 = conductivity_mixture(s_lung_2,s_fluid_2,alpha);
        % Apply left lung first
    lung_elem_indx = fmdl.mat_idx{1,2};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_mix_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_mix_2;
        end
    end
    
    % Apply right lung second
    lung_elem_indx = fmdl.mat_idx{1,3};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
elseif strcmp(pathology, 'ed_right')
    s_mix_1 = conductivity_mixture(s_lung_1,s_fluid_1,alpha);
    s_mix_2 = conductivity_mixture(s_lung_2,s_fluid_2,alpha);
        % Apply left lung first
    lung_elem_indx = fmdl.mat_idx{1,2};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
    % Apply right lung second
    lung_elem_indx = fmdl.mat_idx{1,3};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_mix_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_mix_2;
        end
    end
    
elseif strcmp(pathology, 'ed_low_right')
    s_mix_1 = conductivity_mixture(s_lung_1,s_fluid_1,alpha);
    s_mix_2 = conductivity_mixture(s_lung_2,s_fluid_2,alpha);
    fn_edema = @(x,y,z) ((x-0.45).^2+(y+0.5).^2 < 0.35^2);
    
    % Apply left lung first
    lung_elem_indx = fmdl.mat_idx{1,2};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
    % Apply right lung second
    lung_elem_select = zeros(size(img_sim1.elem_data));
    lung_elem_select(fmdl.mat_idx{1,3}) = true;
    effusion_elem_select = elem_select(img_sim1.fwd_model, fn_edema);
    
    for i=1:length(lung_elem_select) %apply lung area
        if (lung_elem_select(i) == 1) && (img_sim1.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) == 0)
            img_sim1.elem_data(i) = s_lung_1;
        elseif (lung_elem_select(i) == 1) && (img_sim1.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) ~= 0)
            img_sim1.elem_data(i) = s_mix_1*effusion_elem_select(i);
        end
        if (lung_elem_select(i) == 1) && (img_sim2.elem_data(i) ~= (s_heart_2)) && (effusion_elem_select(i) == 0)
            img_sim2.elem_data(i) = s_lung_2;
        elseif (lung_elem_select(i) == 1) && (img_sim2.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) ~= 0)
            img_sim2.elem_data(i) = s_mix_2*effusion_elem_select(i);
        end
    end
    
elseif strcmp(pathology, 'ed_low_left')
    s_mix_1 = conductivity_mixture(s_lung_1,s_fluid_1,alpha);
    s_mix_2 = conductivity_mixture(s_lung_2,s_fluid_2,alpha);
    fn_edema = @(x,y,z) ((x+0.4).^2+(y+0.4).^2 < 0.35^2);
    
    % Apply left lung first
    lung_elem_select = zeros(size(img_sim1.elem_data));
    lung_elem_select(fmdl.mat_idx{1,2}) = true;
    effusion_elem_select = elem_select(img_sim1.fwd_model, fn_edema);
    
    for i=1:length(lung_elem_select) %apply lung area
        if (lung_elem_select(i) == 1) && (img_sim1.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) == 0)
            img_sim1.elem_data(i) = s_lung_1;
        elseif (lung_elem_select(i) == 1) && (img_sim1.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) ~= 0)
            img_sim1.elem_data(i) = s_mix_1*effusion_elem_select(i);
        end
        if (lung_elem_select(i) == 1) && (img_sim2.elem_data(i) ~= (s_heart_2)) && (effusion_elem_select(i) == 0)
            img_sim2.elem_data(i) = s_lung_2;
        elseif (lung_elem_select(i) == 1) && (img_sim2.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) ~= 0)
            img_sim2.elem_data(i) = s_mix_2*effusion_elem_select(i);
        end
    end
    
    % Apply right lung second
    lung_elem_indx = fmdl.mat_idx{1,3};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
elseif strcmp(pathology, 'pleu_eff_right') %fertig
    
    fn_effusion = @(x,y,z) ((x-0.55).^2+(y+0.35).^2 < 0.3^2) | ((x-0.8).^2+(y+0.2).^2 < 0.2^2) | ((x-0.4).^2+(y+0.6).^2 < 0.4^2);
    
    % Apply left lung first
    lung_elem_indx = fmdl.mat_idx{1,2};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
    % Apply right lung second
    lung_elem_select = zeros(size(img_sim1.elem_data));
    lung_elem_select(fmdl.mat_idx{1,3}) = true;
    effusion_elem_select = elem_select(img_sim1.fwd_model, fn_effusion);
    
    for i=1:length(lung_elem_select) %apply lung area
        if (lung_elem_select(i) == 1) && (img_sim1.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) == 0)
            img_sim1.elem_data(i) = s_lung_1;
        elseif (lung_elem_select(i) == 1) && (img_sim1.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) ~= 0)
            img_sim1.elem_data(i) = s_blood_1*effusion_elem_select(i);
        end
        if (lung_elem_select(i) == 1) && (img_sim2.elem_data(i) ~= (s_heart_2)) && (effusion_elem_select(i) == 0)
            img_sim2.elem_data(i) = s_lung_2;
        elseif (lung_elem_select(i) == 1) && (img_sim2.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) ~= 0)
            img_sim2.elem_data(i) = s_blood_2*effusion_elem_select(i);
        end
    end
    
elseif strcmp(pathology, 'pleu_eff_left')
    
    fn_effusion = @(x,y,z) ((x+0.6).^2+(y+0.3).^2 < 0.3^2) | ((x+0.7).^2+(y+0.15).^2 < 0.2^2) | ((x+0.3).^2+(y+0.5).^2 < 0.4^2);
    
    % Apply left lung first
    lung_elem_select = zeros(size(img_sim1.elem_data));
    lung_elem_select(fmdl.mat_idx{1,2}) = true;
    effusion_elem_select = elem_select(img_sim1.fwd_model, fn_effusion);
    
    for i=1:length(lung_elem_select) %apply lung area
        if (lung_elem_select(i) == 1) && (img_sim1.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) == 0)
            img_sim1.elem_data(i) = s_lung_1;
        elseif (lung_elem_select(i) == 1) && (img_sim1.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) ~= 0)
            img_sim1.elem_data(i) = s_blood_1*effusion_elem_select(i);
        end
        if (lung_elem_select(i) == 1) && (img_sim2.elem_data(i) ~= (s_heart_2)) && (effusion_elem_select(i) == 0)
            img_sim2.elem_data(i) = s_lung_2;
        elseif (lung_elem_select(i) == 1) && (img_sim2.elem_data(i) ~= (s_heart_1)) && (effusion_elem_select(i) ~= 0)
            img_sim2.elem_data(i) = s_blood_2*effusion_elem_select(i);
        end
    end
    
    % Apply right lung second
    lung_elem_indx = fmdl.mat_idx{1,3};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
elseif strcmp(pathology, 'pneu_left')
    lung_elem_indx = fmdl.mat_idx{1,2};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end

    llungb_elem_indx = fmdl.mat_idx{1,3};
    llung_elem_indx = fmdl.mat_idx{1,4};
    for i=1:length(llungb_elem_indx) %apply lung area
        if img_sim1.elem_data(llungb_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(llungb_elem_indx(i)) = 0;
        end
        if img_sim2.elem_data(llungb_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(llungb_elem_indx(i)) = 0;
        end
    end
    img_sim1.elem_data(llung_elem_indx) = s_lung_1;
    img_sim2.elem_data(llung_elem_indx) = s_lung_2;

    
else
    
    % Apply left lung first
    lung_elem_indx = fmdl.mat_idx{1,2};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
    
    % Apply right lung second
    lung_elem_indx = fmdl.mat_idx{1,3};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
        end
    end
        
end

% Apply the spine independent of the pathology type
img_sim1.elem_data(spine_indx) = s_spine_1;
img_sim2.elem_data(spine_indx) = s_spine_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% create the difference image
img_diff.elem_data = img_sim2.elem_data-img_sim1.elem_data;


end