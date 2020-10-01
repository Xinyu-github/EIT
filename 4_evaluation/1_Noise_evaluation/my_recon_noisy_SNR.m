function [res_SNR,rms_signal,rms_noise,res_img_noise,img_solve] = my_recon_noisy_SNR(hp,beta,img_sim1,img_sim2,fmdl,SNR_in_dB,img_solve_0)

v1 = fwd_solve(img_sim1); % The voltages at freq_1 are simulated
v2 = fwd_solve(img_sim2); % The voltages at freq_2 are simulated

v1.meas = awgn(v1.meas,SNR_in_dB,'measured');
v2.meas = awgn(v2.meas,SNR_in_dB,'measured');

% weighted_difference 
vL_w1 = v1; vH_w1 = v2; vL_w2 = v1; vH_w2 = v2;

% 2nd WD edition ******
alpha1 = (dot(v1.meas,v2.meas)/dot(v1.meas,v1.meas));

vL_w2.meas = (alpha1).*((v1.meas));
vH_w2.meas = (v2.meas);

hp_param_all = exp((0.5:0.5:100)/2)/(5*10^19); hp_param = hp_param_all(140:200);
% beta = 5; %ii=48;
fmdl_w2 = fmdl;
% weighted diff at measured input V + calculated V
init_iter = 0.01+1i*0.01; %0.01,0.01

weighted_diff = alpha1;
fmdl_w2    = set_weighted_invprob_properties(fmdl_w2, init_iter, hp, weighted_diff,beta);

img_solve = inv_solve(fmdl_w2, vL_w2, vH_w2);
[res_SNR,rms_signal,rms_noise,res_img_noise] = my_eval_resulting_noise(img_solve_0,img_solve);


