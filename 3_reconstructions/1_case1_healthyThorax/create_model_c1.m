function [fmdl,img_sim1,img_sim2,img_diff,org_cond] = create_model_c0(freq_1,freq_2,pathology,elec_num)
%% create FEM of a thorax depending on the pathology
% possible models are healthy, atelectasis, edema

%% create conductivity dist. images
% create conductivities dist. for 2 different freqs
tissues_1 = getTissues_mod(freq_1);
tissues_2 = getTissues_mod(freq_2); 
% create image included organs with muscle background for simulation
s_bkgnd_1 = round(tissues_1.Muscle_combined,4,'decimal');
s_bkgnd_2 = round(tissues_2.Muscle_combined,4,'decimal'); 
%get conductivites
s_lung_1 = round(tissues_1.LungInflated_combined,4,'decimal'); 
s_heart_1 = round(tissues_1.Heart_combined,4,'decimal'); 
s_spine_1 = round(tissues_1.BoneAussen_combined,4,'decimal');
s_fluid_1 = round(tissues_1.fluid_combined,4,'decimal'); 
s_lung_2 = round(tissues_2.LungInflated_combined,4,'decimal'); 
s_heart_2 = round(tissues_2.Heart_combined,4,'decimal'); 
s_spine_2 = round(tissues_2.BoneAussen_combined,4,'decimal');
s_fluid_2 = round(tissues_1.fluid_combined,4,'decimal'); 

if pathology == 1 %healthy
    fmdl = create_thorax_fem(elec_num,3,0.03,1,0);
    s_lung_right_1=s_lung_1;
    s_lung_right_2=s_lung_2;
elseif pathology == 2 %atelectasis
    fmdl = create_thorax_fem(elec_num,3,0.03,2,0);
elseif pathology == 3 %edema
    fmdl = create_thorax_fem(elec_num,3,0.03,1,0);
    s_lung_right_1=s_fluid_1;
    s_lung_right_2=s_fluid_2;
else
    fmdl = create_thorax_fem(elec_num,3,0.03,4,0);
end

% create image included organs with muscle background for simulation

fmdl.organs.bg1 = s_bkgnd_1; fmdl.organs.bg2 = s_bkgnd_2;
fmdl.organs.heart1 = s_heart_1; fmdl.organs.heart2 = s_heart_2;
fmdl.organs.lung1 = s_lung_1; fmdl.organs.lung2 = s_lung_2;
fmdl.organs.spine1 = s_spine_1; fmdl.organs.spine2 = s_spine_2;

img_sim1 = mk_image( fmdl, s_bkgnd_1 );  
img_sim2 = mk_image( fmdl, s_bkgnd_2 ); 
img_diff = img_sim1;

%% create organ shapes
heart_center = [0.12,0.25]; % assume circle heart, to place position of the heart, 
spine_center = [0.01,-0.47]; % assume ellipse, to place position of the spine
lung_center_right = [0.15,0.1]; % assume cropped circle, to place position of the right lung (the left lung is mirrored)
heart_radius = 0.27;% Radius of the heart (modeled by a circle)
organ_fn = create_simple_organ_shape( heart_center,heart_radius, spine_center, lung_center_right );

% extract conductivity value from database and function to create organ shape
fn_heart = organ_fn.heart;
fn_spine = organ_fn.spine;
%fn_lung = organ_fn.lung;

org_cond.bg1 = s_bkgnd_1; org_cond.bg2 = s_bkgnd_2;
org_cond.heart1 = s_heart_1; org_cond.heart2 = s_heart_2;
org_cond.lung1 = s_lung_1; org_cond.lung2 = s_lung_2;
org_cond.spine1 = s_spine_1; org_cond.spine2 = s_spine_2;


%% apply conductivities on each organ 
img_sim1.elem_data(:) = s_bkgnd_1 ... %apply heart, spine
                +( (s_heart_1)- (s_bkgnd_1) )*elem_select(img_sim1.fwd_model, fn_heart)...
                +( (s_spine_1)- (s_bkgnd_1) )*elem_select(img_sim1.fwd_model, fn_spine);
img_sim2.elem_data(:) = s_bkgnd_2 ... %apply heart, spine
                +( (s_heart_2)- (s_bkgnd_2) )*elem_select(img_sim2.fwd_model, fn_heart)...
                +( (s_spine_2)- (s_bkgnd_2) )*elem_select(img_sim2.fwd_model, fn_spine);
spine_indx = find(round(img_sim1.elem_data(:),4,'decimal')==round(s_spine_1,4,'decimal'));
%%%%%  apply lung area  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
lung_elem_indx = fmdl.mat_idx{1,2};
for i=1:length(lung_elem_indx) %apply lung area
    if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
        img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
    end
    if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
        img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
    end
end
if pathology ~= 2
    lung_elem_indx = fmdl.mat_idx{1,3};
    for i=1:length(lung_elem_indx) %apply lung area
        if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
            img_sim1.elem_data(lung_elem_indx(i)) = s_lung_right_1;
        end
        if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
            img_sim2.elem_data(lung_elem_indx(i)) = s_lung_right_2;
        end
    end
end
img_sim1.elem_data(spine_indx) = s_spine_1;
img_sim2.elem_data(spine_indx) = s_spine_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
img_diff.elem_data = img_sim2.elem_data-img_sim1.elem_data;
end
