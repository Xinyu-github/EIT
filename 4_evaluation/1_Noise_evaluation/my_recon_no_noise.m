function [img_solve_w2_0] = my_recon_no_noise(hp,beta,img_sim1,img_sim2,fmdl)
%% no noise
v1 = fwd_solve(img_sim1); % The voltages at freq_1 are simulated
v2 = fwd_solve(img_sim2); % The voltages at freq_2 are simulated
% weighted_difference 
vL_w1 = v1; vH_w1 = v2; vL_w2 = v1; vH_w2 = v2;

% 2nd WD edition ******
alpha1 = (dot(v1.meas,v2.meas)/dot(v1.meas,v1.meas));

vL_w2.meas = (alpha1).*((v1.meas));
vH_w2.meas = (v2.meas);

% beta = 5; 
fmdl_w2 = fmdl;
% weighted diff at measured input V + calculated V
init_iter = 0.01+1i*0.01; %0.01,0.01

weighted_diff = alpha1;
fmdl_w2    = my_weighted_invprob_properties(fmdl_w2, init_iter, hp, weighted_diff,beta);

img_solve_w2_0= inv_solve(fmdl_w2, vL_w2, vH_w2);





