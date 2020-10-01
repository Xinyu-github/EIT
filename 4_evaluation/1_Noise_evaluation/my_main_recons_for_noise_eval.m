clc,clear;
%% create FEM
fmdl = my_create_thorax_fem(32,0.03);
%% create conductivity dist. images
% case 4: complex model of cross-section of body with organ:heart,left-lung,right-lung,spine
% Select 2 injected frequencies (Hz)
freq_1 = 100e3; freq_2 = 700e3; 
% create conductivities dist. for 2 different freqs
tissues_1 = my_getTissues(freq_1);% tissues_1 = my_create_tissue_conductivity(freq_1); 
tissues_2 = my_getTissues(freq_2);% tissues_2 = my_create_tissue_conductivity(freq_2); 
% create image included organs with muscle background for simulation
s_bkgnd_1 = tissues_1.Muscle_combined; s_bkgnd_2 = tissues_2.Muscle_combined; 
img_sim1 = mk_image( fmdl, s_bkgnd_1 );  
img_sim2 = mk_image( fmdl, s_bkgnd_2 ); 
img_diff = img_sim1; %figure,show_fem(img_diff)

%% create organ shapes
heart_center = [0.62,0.72]; % assume circle heart, to place position of the heart, 
spine_center = [0.55,0.1]; % assume ellipse, to place position of the spine
lung_center_right = [0.15,0.1]; % assume cropped circle, to place position of the right lung (the left lung is mirrored)
heart_radius = 0.3;% Radius of the heart (modeled by a circle)
organ_fn = my_create_simple_organ_shape( heart_center,heart_radius, spine_center);

% extract conductivity value from database and function to create organ shape
fn_heart = organ_fn.heart; fn_spine = organ_fn.spine; 
s_lung_1 = tissues_1.LungInflated_combined; s_heart_1 = tissues_1.Heart_combined; 
s_spine_1 = tissues_1.BoneAussen_combined; 
s_lung_2 = tissues_2.LungInflated_combined; s_heart_2 = tissues_2.Heart_combined; 
s_spine_2 = tissues_2.BoneAussen_combined; 

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
lung_elem_indx = fmdl.mat_idx{1,3};
for i=1:length(lung_elem_indx) %apply lung area
    if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
        img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
    end
    if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
        img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
    end
end
img_sim1.elem_data(spine_indx) = s_spine_1;
img_sim2.elem_data(spine_indx) = s_spine_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
img_diff.elem_data = img_sim2.elem_data-img_sim1.elem_data;
hp_param_all = exp((0.5:0.5:100)/2)/(5*10^19); hp_param = hp_param_all(140:200);

%%
%hp=8.5117; B_mixed=5; %these hp and B can be modified!!!!!!!!!!!
hp=31.2321; B_mixed=0.5;
% Reconstruction without noise
[img_solve_0] = my_recon_no_noise(hp,B_mixed,img_sim1,img_sim2,fmdl);
save('img_solve_hp31_snr0','img_solve_0');
%with 0.5beta

for ind = 1:length(noise)
    SNR_in_dB = noise(ind);
    [res_SNR,rms_signal,rms_noise,res_img_noise,img_solve] = my_recon_noisy_SNR(hp,B_mixed,img_sim1,img_sim2,fmdl,SNR_in_dB,img_solve_0);
    RMS_s(ind) = rms_signal;
    RMS_n(ind) = rms_noise;
    res_dB_n(ind) = res_SNR;
    img_noise.img(ind) = res_img_noise;
    img_noisy_res.img(ind) = img_solve;
end
save('RMS_s_hp31_6','RMS_s');
save('RMS_n_hp31_6','RMS_n');
save('res_dB_n_hp31_6','res_dB_n');
save('img_noise_hp31_6','img_noise');
save('img_noisy_res_hp31_6','img_noisy_res');


% Reconstruction with noise
% Due to random noise, 5 datasets were simulated
noise = [0:5:100];
for ind = 1:length(noise)
    SNR_in_dB = noise(ind);
    [res_SNR,rms_signal,rms_noise,res_img_noise,img_solve] = my_recon_noisy_SNR(hp,B_mixed,img_sim1,img_sim2,fmdl,SNR_in_dB,img_solve_0);
    RMS_s(ind) = rms_signal;
    RMS_n(ind) = rms_noise;
    res_dB_n(ind) = res_SNR;
    img_noise.img(ind) = res_img_noise;
    img_noisy_res.img(ind) = img_solve;
end
save('RMS_s_hp851_1','RMS_s');
save('RMS_n_hp851_1','RMS_n');
save('res_dB_n_hp851_1','res_dB_n');
save('img_noise_hp851_1','img_noise');
save('img_noisy_res_hp851_1','img_noisy_res');

for ind = 1:length(noise)
    SNR_in_dB = noise(ind);
    [res_SNR,rms_signal,rms_noise,res_img_noise,img_solve] = my_recon_noisy_SNR(hp,B_mixed,img_sim1,img_sim2,fmdl,SNR_in_dB,img_solve_0);
    RMS_s(ind) = rms_signal;
    RMS_n(ind) = rms_noise;
    res_dB_n(ind) = res_SNR;
    img_noise.img(ind) = res_img_noise;
    img_noisy_res.img(ind) = img_solve;  
end
save('RMS_s_hp851_2','RMS_s');
save('RMS_n_hp851_2','RMS_n');
save('res_dB_n_hp851_2','res_dB_n');
save('img_noise_hp851_2','img_noise');
save('img_noisy_res_hp851_2','img_noisy_res');

for ind = 1:length(noise)
    SNR_in_dB = noise(ind);
    [res_SNR,rms_signal,rms_noise,res_img_noise,img_solve] = my_recon_noisy_SNR(hp,B_mixed,img_sim1,img_sim2,fmdl,SNR_in_dB,img_solve_0);
    RMS_s(ind) = rms_signal;
    RMS_n(ind) = rms_noise;
    res_dB_n(ind) = res_SNR;
    img_noise.img(ind) = res_img_noise;
    img_noisy_res.img(ind) = img_solve;  
end
save('RMS_s_hp851_3','RMS_s');
save('RMS_n_hp851_3','RMS_n');
save('res_dB_n_hp851_3','res_dB_n');
save('img_noise_hp851_3','img_noise');
save('img_noisy_res_hp851_3','img_noisy_res');

for ind = 1:length(noise)
    SNR_in_dB = noise(ind);
    [res_SNR,rms_signal,rms_noise,res_img_noise,img_solve] = my_recon_noisy_SNR(hp,B_mixed,img_sim1,img_sim2,fmdl,SNR_in_dB,img_solve_0);
    RMS_s(ind) = rms_signal;
    RMS_n(ind) = rms_noise;
    res_dB_n(ind) = res_SNR;
    img_noise.img(ind) = res_img_noise;
    img_noisy_res.img(ind) = img_solve;   
end
save('RMS_s_hp851_4','RMS_s');
save('RMS_n_hp851_4','RMS_n');
save('res_dB_n_hp851_4','res_dB_n');
save('img_noise_hp851_4','img_noise');
save('img_noisy_res_hp851_4','img_noisy_res');

for ind = 1:length(noise)
    SNR_in_dB = noise(ind);
    [res_SNR,rms_signal,rms_noise,res_img_noise,img_solve] = my_recon_noisy_SNR(hp,B_mixed,img_sim1,img_sim2,fmdl,SNR_in_dB,img_solve_0);
    RMS_s(ind) = rms_signal;
    RMS_n(ind) = rms_noise;
    res_dB_n(ind) = res_SNR;
    img_noise.img(ind) = res_img_noise;
    img_noisy_res.img(ind) = img_solve;
end
save('RMS_s_hp851_5','RMS_s');
save('RMS_n_hp851_5','RMS_n');
save('res_dB_n_hp851_5','res_dB_n');
save('img_noise_hp851_5','img_noise');
save('img_noisy_res_hp851_5','img_noisy_res');