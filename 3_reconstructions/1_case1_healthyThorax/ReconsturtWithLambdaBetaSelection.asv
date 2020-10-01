%% 1) Create input model of conductivity distribution 
%%% -Model with conductivity assigned for higher and lower frequencies,
%%% img_sim2 and img_sim1, respectively.
%%% -img_diff is the conductivity difference between img_sim2 and img_sim1
%%% -fmdl is forward model
%%% -org collects conductivity for all organs
nSkip = 0;
freq_1 = 100e3; freq_2 = 700e3;
[fmdl,img_sim1,img_sim2,img_diff] = create_model_pat(freq_1, freq_2, 'healthy', 32, 0.03, nSkip);
%plot_fdEITmodel(img_diff, 'model_c1', false);

%%plot model for publication
%figure; show_fem(fmdl, [0 1]); axis off; savefig('mesh_thorax_axis');

%% 2) Generate voltage measurement
%%% Calculate the measured voltages via solving forward problem
v1 = fwd_solve(img_sim1); % The voltages at freq_1 are simulated
v2 = fwd_solve(img_sim2); % The voltages at freq_2 are simulated
%%% Add input noise
% SNR_in_dB = 30;
% v1.meas = awgn(v1.meas,SNR_in_dB,'measured');
% v2.meas = awgn(v2.meas,SNR_in_dB,'measured');

init_iter = 0.5+1i*0.02; % initial conductivity
fmdl_recon = create_thorax_fem_simon(32, nSkip, 0.05, 0, 0, 0);
%figure,show_fem(mk_image(fmdl_recon,1),[0,1]);
imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, 0,0,true,0); %initialize inverse model
init_iter = init_optimization(imdl,v1,v2,'absolute'); %optimize initial guess 
%% calculate electrode distance
% dist=0;
% for jj=1:length(fmdl.electrode)
%     dist = dist + norm(mean(fmdl_recon.nodes(fmdl_recon.electrode(jj).nodes(:),:))-mean(fmdl.nodes(fmdl.electrode(jj).nodes(:),:)));
% end
% meanDist = dist/length(fmdl.electrode) %0.0664 for angular even spacing; own spacing: 0.0431
% electrode_error = meanDist/(max(fmdl_recon.nodes(:,2))-min(fmdl_recon.nodes(:,2))) %0.0499 for angular even spacing; own spacing: 0.0324

%% start the reconstruction and store under /results
hp_param=logspace(-6,-1,40);
mode= 1; %choosing different objective function
% Reconstruction in weighted differential form
for ii = 1:length(hp_param)
    imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, hp_param(ii),1,0,mode);
    img_mode1(ii)= inv_solve(imdl, v1, v2);

    disp('o-----------------------------------------------------------------------------------o')
    disp(strcat('Finished reconstruction for hyperparameter: '));
    disp(strcat(num2str(ii),'/',num2str(length(hp_param)),' -> hp=',num2str(hp_param(ii))));
    disp('o-----------------------------------------------------------------------------------o')
end
save(strcat('./results/recon_c1_mode1_allHP.mat'),'img_mode1','hp_param');

% % Reconstruction in seperated absolute form
mode= 2;
for ii = 1:length(hp_param)
    imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, hp_param(ii),1,0,mode);
    img_mode2(ii)= inv_solve(imdl, v1, v2);

    disp('o-----------------------------------------------------------------------------------o')
    disp(strcat('Finished reconstruction for hyperparameter: '));
    disp(strcat(num2str(ii),'/',num2str(length(hp_param)),' -> hp=',num2str(hp_param(ii))));
    disp('o-----------------------------------------------------------------------------------o')
end
save(strcat('./results/recon_c1_mode2_allHP.mat'),'img_mode2','hp_param');

%% Selection of suitable hyperparameters for both reconstruction

% L-curve
[res1,reg1]=generate_lcurve(img_mode1);
[res2,reg2]=generate_lcurve(img_mode2);

% % SNR metrics
% mode = 1;
% imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, 0,1,0,mode);
% imdl.hyperparameter.n_targets = 250;
% for i = 1: length(hp_param)
%     SNR1(i) = my_calc_image_SNR(imdl,hp_param(i),false,1);
% end
% idx_hpsnr1 = find(SNR1==max(SNR1));
% hp_snr1 = hp_param(idx_hpsnr1);
% 
% mode = 2;
% imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, 1,1,0,mode);
% imdl.hyperparameter.n_targets = 250;
% for i = 1: length(hp_param)
%     SNR2(i) = my_calc_image_SNR(imdl,hp_param(i),false,1);
% end
% idx_hpsnr2 = find(SNR2==max(SNR2));
hp_snr2 = hp_param(idx_hpsnr2);

%% Beta Selection
% ImgDiff = img_mode1(idx_hpsnr1); % best reconstruction obtained by differential norm
% ImgAbs = img_mode2(idx_hpsnr2); % best reconstruction obtained by absolute norm
% Img = beta_selection(ImgDiff,ImgAbs,0); % Padded best result with the best beta
% plot_fdEITmodel(Img,'case1',0);



