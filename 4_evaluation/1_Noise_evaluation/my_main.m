clc,clear;
%% create FEM
n_elec = 32; max_gs = 0.03; %number of electrodes and maximum grid size
fmdl = my_create_thorax_fem(n_elec,max_gs); %create forward model
%% create conductivity dist. images
% Cross-section model of body with organ:heart,left-lung,right-lung,spine
% Select 2 frequencies of applied current (Hz)
freq_1 = 100e3; freq_2 = 700e3; 
% create conductivities dist. for 2 different freqs
tissues_1 = my_getTissues(freq_1);% get tissues conductivities  
tissues_2 = my_getTissues(freq_2);% get tissues conductivities  
% create image included organs with muscle background for simulation
s_bkgnd_1 = tissues_1.Muscle_combined; s_bkgnd_2 = tissues_2.Muscle_combined; %background conductivity
img_sim1 = mk_image( fmdl, s_bkgnd_1 );  %at f1
img_sim2 = mk_image( fmdl, s_bkgnd_2 );  %at f2
img_diff = img_sim1; 

%% create organ shapes
heart_center = [0.62,0.72]; % assume circle heart, to place position of the heart, 
spine_center = [0.55,0.1]; % assume ellipse, to place position of the spine
heart_radius = 0.3;% Radius of the heart (modeled by a circle)
organ_fn = my_create_simple_organ_shape( heart_center, heart_radius, spine_center);

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
lung_elem_indx = fmdl.mat_idx{1,1};
for i=1:length(lung_elem_indx) %apply lung area
    if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
        img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
    end
    if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
        img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
    end
end
% lung_elem_indx = fmdl.mat_idx{1,3};
% for i=1:length(lung_elem_indx) %apply lung area
%     if img_sim1.elem_data(lung_elem_indx(i)) ~= (s_heart_1)
%         img_sim1.elem_data(lung_elem_indx(i)) = s_lung_1;
%     end
%     if img_sim2.elem_data(lung_elem_indx(i)) ~= (s_heart_2)
%         img_sim2.elem_data(lung_elem_indx(i)) = s_lung_2;
%     end
% end
img_sim1.elem_data(spine_indx) = s_spine_1;
img_sim2.elem_data(spine_indx) = s_spine_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
img_diff.elem_data = img_sim2.elem_data-img_sim1.elem_data;


%% forward solve
SNR_in_dB = 40; 
v1 = fwd_solve(img_sim1); % The voltages at freq_1 are simulated
v2 = fwd_solve(img_sim2); % The voltages at freq_2 are simulated

v1.meas = awgn(v1.meas,SNR_in_dB,'measured');
v2.meas = awgn(v2.meas,SNR_in_dB,'measured');

%% weighted_difference 
vL_w2 = v1; vH_w2 = v2;

% Weighted different algorithm ******
alpha1 = (dot(v1.meas,v2.meas)/dot(v1.meas,v1.meas));

vL_w2.meas = (alpha1).*((v1.meas));
vH_w2.meas = (v2.meas);

vL_w2.elem_data = img_sim1.elem_data;
vH_w2.elem_data = img_sim2.elem_data;

%% Reconstruction (solve inverse problem)
hp_param_all = exp((0.5:0.5:100)/2)/(5*10^19); hp_param = hp_param_all(140:200);
beta = 5; %default beta for absolute-mixed minimisation algorithm
hp_param = 6.62; %already selected hp_param by L-curve for this model

% weighted diff at measured input V + calculated V
init_iter = 0.01+1i*0.01; %0.01,0.01

weighted_diff = alpha1;
fmdl_w2    = my_weighted_invprob_properties(fmdl, init_iter, hp_param, weighted_diff,beta);

img_solve_w2= inv_solve(fmdl_w2, vL_w2, vH_w2); 

cm = load("cm_hot.mat"); cm=cm.cm;
figure, show_phase_contour_slice(img_diff,64,50,cm);
figure, show_phase_contour_slice(img_solve_w2,64,50,cm);

%% load saved files of reconstructed and reference images
% save('img_diff','img_diff'), save('img_solve_w2','img_solve_w2')
% img_diff = load('img_diff'); img_solve_w2 = load('img_solve_w2');
% img_diff = img_diff.img_diff; img_solve_w2 = img_solve_w2.img_solve_w2;
% figure, show_phase_contour_slice(img_diff,64,50,cm);
% figure, show_phase_contour_slice(img_solve_w2,64,50,cm);
