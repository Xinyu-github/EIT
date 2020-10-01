%% 1) Create input model of conductivity distribution 
%%% -Model with conductivity assigned for higher and lower frequencies,
%%% img_sim2 and img_sim1, respectively.
%%% -img_diff is the conductivity difference between img_sim2 and img_sim1
%%% -fmdl is forward model
%%% -org collects conductivity for all organs
nSkip = 2;
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

init_iter = 0.1+1i*0.01; % initial conductivity
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
hp_param=logspace(-5,-1,30);
beta=0.*ones(1,length(hp_param));
mode= 1;
img_mode1=reconstruct_multipleHP(fmdl_recon,hp_param,beta,true,init_iter, v1, v2, mode,'abc');

hp_param=logspace(-5,-1,30);
beta=1.*ones(1,length(hp_param));
mode= 2;
img_mode=reconstruct_multipleHP(fmdl_recon,hp_param,beta,true,init_iter, v1, v2, mode,'c1_mode2');
% only absolute term
% 
% % Select variant hyper-parameters and different absolute-mixed parameters
% hp_param=logspace(-7,-1,30);
% beta=[0 0.1 0.25 0.5 1 2.5 5 10]'.*ones(8,length(hp_param));
% reconstruct_multipleHP(fmdl_recon,hp_param,beta,true,init_iter, v1, v2, 'c1_B');
% 
% % no beta no weighted
% hp_param=logspace(-4,-1,40);
% beta=zeros(1,length(hp_param));
% reconstruct_multipleHP(fmdl_recon,hp_param,beta,false,init_iter, v1, v2, 'c1_onlyGN_B');
% 
% % beta is scaled with hp
% hp_param=logspace(-4,-1,40);
% beta=[0.5,1,2,4,6,8]'.*hp_param;
% reconstruct_multipleHP(fmdl_recon,hp_param,beta,true,init_iter, v1, v2, '_c1_scaling_B');
% 
% % beta is double hp
% hp_param=logspace(-4,-1,40);
% beta=2*hp_param;
% reconstruct_multipleHP(fmdl_recon,hp_param,beta,true,init_iter, v1, v2, '_c1_scaling2_increasedResolution_B');
% 
% % beta is 0
% hp_param=logspace(-4,-1,20);
% beta=0*hp_param;
% 
% reconstruct_multipleHP(fmdl_recon,hp_param,beta,true,init_iter, v1, v2, '_c1_3iter');
% 

img = beta_selection(img_mode1(31),img_mode(23));
