nSkip = 0;

% Specify the simulation frequencies
freq_1 = 100e3; 
freq_2 = 700e3;

beta= 0.1;
init_iter = 0.1+1i*0.01; 
%% Inititalize the cell array to store the models and their respective names
cases = cell(1,11);

[cases{1}.fmdl,cases{1}.img_sim1,cases{1}.img_sim2,cases{1}.img_diff] = create_model_pat(freq_1, freq_2, 'healthy', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{1}.img_diff, 'model_healthy', 0,'t');
cases{1}.title = 'healthy';

[cases{2}.fmdl,cases{2}.img_sim1,cases{2}.img_sim2,cases{2}.img_diff] = create_model_pat(freq_1, freq_2, 'at_right', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{2}.img_diff, 'model_at_right', 0,'t');
cases{2}.title = 'at_right';

[cases{3}.fmdl,cases{3}.img_sim1,cases{3}.img_sim2,cases{3}.img_diff] = create_model_pat(freq_1, freq_2, 'at_left', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{3}.img_diff, 'model_at_left', 0,'t');
cases{3}.title = 'at_left';

[cases{4}.fmdl,cases{4}.img_sim1,cases{4}.img_sim2,cases{4}.img_diff] = create_model_pat(freq_1, freq_2, 'ed_right', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{4}.img_diff, 'model_ed_right', 0,'t');
cases{4}.title = 'ed_right';

[cases{5}.fmdl,cases{5}.img_sim1,cases{5}.img_sim2,cases{5}.img_diff] = create_model_pat(freq_1, freq_2, 'ed_left', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{5}.img_diff, 'model_ed_left', 0,'t');
cases{5}.title = 'ed_left';

[cases{6}.fmdl,cases{6}.img_sim1,cases{6}.img_sim2,cases{6}.img_diff] = create_model_pat(freq_1, freq_2, 'peri_eff', 32, 0.03, nSkip);
% plot_fdEITmodel(cases{6}.img_diff, 'model_peri_eff', 0,'t');
cases{6}.title = 'peri_eff';

[cases{7}.fmdl,cases{7}.img_sim1,cases{7}.img_sim2,cases{7}.img_diff] = create_model_pat(freq_1, freq_2, 'ed_low_right', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{7}.img_diff, 'model_ed_low_right', 0,'t');
cases{7}.title = 'ed_low_right';

[cases{8}.fmdl,cases{8}.img_sim1,cases{8}.img_sim2,cases{8}.img_diff] = create_model_pat(freq_1, freq_2, 'ed_low_left', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{8}.img_diff, 'model_ed_low_left', 0,'t');
cases{8}.title = 'ed_low_left';

[cases{9}.fmdl,cases{9}.img_sim1,cases{9}.img_sim2,cases{9}.img_diff] = create_model_pat(freq_1, freq_2, 'pleu_eff_right', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{9}.img_diff, 'model_pleu_eff_right', 0,'t');
cases{9}.title = 'pleu_eff_right';

[cases{10}.fmdl,cases{10}.img_sim1,cases{10}.img_sim2,cases{10}.img_diff] = create_model_pat(freq_1, freq_2, 'pleu_eff_left', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{10}.img_diff, 'model_pleu_eff_left', 0,'t');
cases{10}.title = 'pleu_eff_left';

[cases{11}.fmdl,cases{11}.img_sim1,cases{11}.img_sim2,cases{11}.img_diff] = create_model_pat(freq_1, freq_2, 'pneu_left', 32, 0.03, nSkip);
%plot_fdEITmodel(cases{11}.img_diff, 'model_pneu_left', 0,'t');
cases{11}.title = 'pneu_left';
%% Calculate the measured voltages via solving forward problem
for i = 1:length(cases)
cases{i}.v1 = fwd_solve(cases{i}.img_sim1); % The voltages at freq_1 are simulated
cases{i}.v2 = fwd_solve(cases{i}.img_sim2); % The voltages at freq_2 are simulated
end

fmdl_recon = create_thorax_fem_simon(32, nSkip, 0.05, 0, 0, 0);
%% Selection of hyperparameter


%%% SNR metric
imdl    = mk_hpselection_imdl(cases{1}.fmdl,fmdl_recon, 0, beta,'snr');
% find hyperparameter which can obtain a minimal snr 
hp_snr = fminbnd(@(x)my_calc_image_SNR(imdl,x,true),0,0.01);


%%% BestRes
imdl    = mk_hpselection_imdl(cases{1}.fmdl,fmdl_recon, 0, beta,'bestres');
hp_br = fminbnd(@(x)calc_image_br(imdl,x),0,1);

% % % Lcurve
v1 = cases{4}.v1.meas;
v2 = cases{4}.v2.meas;



% Initialize the range of hyperparameter
Hp = [5e-3:1e-4:9e-3];
% Reconstrut for each hp one image

k = 1; % specify the case that to be reconstructed
imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, 0,true,beta); %initialize inverse model
optinit = init_optimization(imdl,cases{k}.v1,cases{k}.v2,'absolute'); %optimize initial guess 
clear img_solve
for ii = 1:length(Hp)
    imdl    = set_weighted_invprob_properties(fmdl_recon, optinit, Hp(ii), true,beta);
    img_solve(ii)= inv_solve(imdl, cases{k}.v1.meas, cases{k}.v2.meas);
end

%%% L-curve

[res,reg] = generate_lcurve(img_step);
