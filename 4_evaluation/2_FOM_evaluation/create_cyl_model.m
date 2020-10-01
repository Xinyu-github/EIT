function [v1, v2, v1w, v2w, alpha1, imdl, img_diff] = create_cyl_model(cyl_center)
%%%%% input: cyl_center
%%%%% output: v1w,v2w
%imdl = load('fmdl_cyl.mat'); imdl = imdl.fmdl;
imdl = mk_common_model('h2C');
%show_fem(mk_image(imdl,1));
freq_1 = 100e3; freq_2 = 700e3; freq_0 = 1e2;
% create conductivities dist. for 2 different freqs
tissues_1 = getTissues_mod(freq_1); tissues_2 = getTissues_mod(freq_2);
% create image included organs with muscle background for simulation
s_bkgnd_1 = round(tissues_1.Muscle_combined,4,'decimal'); 
s_bkgnd_2 = round(tissues_2.Muscle_combined,4,'decimal'); 
s_lung_1 = round(tissues_1.LungInflated_combined,4,'decimal'); 
s_lung_2 = round(tissues_2.LungInflated_combined,4,'decimal'); 
img_sim1 = mk_image( imdl, s_bkgnd_1 );  
img_sim2 = mk_image( imdl, s_bkgnd_2 ); 
img_diff = img_sim1; %figure,show_fem(img_diff)
%% create cyl shapes with Inflated lung cond and bg with muscle cond
radius_max=max(sqrt(img_diff.fwd_model.nodes(img_diff.fwd_model.boundary(:,2),1).^2+img_diff.fwd_model.nodes(img_diff.fwd_model.boundary(:,2),2).^2));
cyl_radius = 0.05*2*radius_max; % Radius of cyl
fn_cyl = @(x,y,z) (x-cyl_center(1)).^2+(y-cyl_center(2)).^2 < cyl_radius^2;

%% apply conductivities on each organ 
img_sim1.elem_data(:) = s_bkgnd_1 ... %apply heart, spine
                +( (s_lung_1)- (s_bkgnd_1) )*elem_select(img_sim1.fwd_model, fn_cyl);
% for homogenous background over frequency
%s_bkgnd_2 = s_bkgnd_1;
img_sim2.elem_data(:) = s_bkgnd_2 ... %apply heart, spine
                +( (s_lung_2)- (s_bkgnd_2) )*elem_select(img_sim2.fwd_model, fn_cyl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
img_diff.elem_data = img_sim2.elem_data-img_sim1.elem_data;

%relative_change=(s_lung_2)-(s_lung_1)/(s_bkgnd_2-s_bkgnd_1);

% figure, show_fem(img_diff),title('X_t');

%% forward solve
%SNR_in_dB = 40; 
v1 = fwd_solve(img_sim1); % The voltages at freq_1 are simulated
v2 = fwd_solve(img_sim2); % The voltages at freq_2 are simulated

%v1.meas = awgn(v1.meas,SNR_in_dB,'measured');
%v2.meas = awgn(v2.meas,SNR_in_dB,'measured');

%% weighted_difference 
v1w = v1; v2w = v2;

% 2nd WD edition ******
alpha1 = (dot(v1.meas,v2.meas)/dot(v1.meas,v1.meas));

v1w.meas = (alpha1).*((v1.meas));
v2w.meas = (v2.meas);

