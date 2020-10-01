freq_1 = 100e3;
freq_2 = 700e3;

%% Generate model
% 2D Model
imdl_2d= mk_common_model('d2d2c',32);
img_sim1= mk_image(imdl_2d.fwd_model,1);
img_sim2 = img_sim1;
img_diff = img_sim1;

tissues_1 = getTissues_mod(freq_1);
tissues_2 = getTissues_mod(freq_2); 

% Get muscle conductivity as background at both frequencies
s_bkgnd_1 = round(tissues_1.Muscle_combined,4,'decimal'); 
s_bkgnd_2 = round(tissues_2.Muscle_combined,4,'decimal'); 

% Get conductivites at frequency 1
s_lung_1 = round(tissues_1.LungInflated_combined,4,'decimal'); 

% Get conductivitied at frequency 2
s_lung_2 = round(tissues_2.LungInflated_combined,4,'decimal'); 

select_fcn = inline('(x+0.8).^2+(y-0).^2<0.15^2','x','y','z');
img_sim1.elem_data = s_bkgnd_1 + (s_lung_1-s_bkgnd_1)*elem_select(img_sim1.fwd_model, select_fcn);

img_sim2.elem_data = s_bkgnd_2 + (s_lung_2-s_bkgnd_2)*elem_select(img_sim2.fwd_model, select_fcn);

img_diff.elem_data = img_sim2.elem_data - img_sim1.elem_data;

v1 = fwd_solve(img_sim1);
v2 = fwd_solve(img_sim2);
%% Selection of hyperparameter
% beta= 0;
init_iter = 0.5+0.02i; 
imdl = mk_common_model('c2d2c',32);
fmdl_recon = imdl.fwd_model;
imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, 0,0,true,1); %initialize inverse model
init_iter = init_optimization(imdl,v1,v2,'absolute'); %optimize initial guess 
% %%% SNR metric
% imdl    = my_weighted_invprob_properties(fmdl.fwd_model, init_iter, 0, true,beta);
% imdl.hyperparameter.n_targets = 500;
% % Specify noise
% NOISE_AMPLITUDE = 1E-3;
% imdl.noise1 = NOISE_AMPLITUDE * randn(32*32, 1000);
% imdl.noise2 = NOISE_AMPLITUDE * randn(32*32, 1000);

% find hyperparameter which can obtain a minimal snr 
% hp_snr = fminbnd(@(x)my_calc_image_SNR(imdl,x,true),0,1);
%% 
% 
% hp = logspace(-6,0,20);
% beta = 0;
% BR = zeros(length(beta),length(hp));
% for i = 1:length(beta)
%     for j = 1:length(hp)
%         imdl    = mk_hpselection_imdl(imdl_2d,fmdl_recon, hp(j), beta(i),0,'bestres');
%         BR(i,j) = calc_image_br(imdl,hp(j));
%     end
% end
%%% L-curve
Hp1 = logspace(-6,-1,20);
Hp2 = logspace(-6,-1,20);

for ii = 1:length(Hp1)
    imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, Hp1(ii),1,...
        1,1);
    img1(ii)= inv_solve(imdl, v1, v2);
    disp('o-----------------------------------------------------------------------------------o')
   
    disp(strcat(num2str(ii),'/',num2str(length(Hp1)),' -> hp=',num2str(Hp1(ii))));
    disp('o-----------------------------------------------------------------------------------o')
end
for ii = 1:length(Hp2)
    imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, Hp2(ii),1,...
        1,2);
    img2(ii)= inv_solve(imdl, v1, v2);

    disp('o-----------------------------------------------------------------------------------o')
   
    disp(strcat(num2str(ii),'/',num2str(length(Hp2)),' -> hp=',num2str(Hp2(ii))));
    disp('o-----------------------------------------------------------------------------------o')
end
% [res,reg] = generate_lcurve(imgr1);
img_beta = cell(1,length(beta));
for jj =1:length(beta)
    for ii = 1:length(Hp1)
    imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, Hp1(ii),beta(jj),...
        1,1);
    img(ii)= inv_solve(imdl, v1, v2);
    disp('o-----------------------------------------------------------------------------------o')
   
    disp(strcat(num2str(ii),'/',num2str(length(Hp1)),' -> hp=',num2str(Hp1(ii))));
    disp(strcat(num2str(ii),'/',num2str(length(beta)),' -> beta=',num2str(beta(jj))));
    disp('o-----------------------------------------------------------------------------------o');
    img_beta{jj} = img;
    end
end