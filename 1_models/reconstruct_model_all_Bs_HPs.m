% This script can be used to generate 10 different thorax models with
% different lung pathologies. These models are then used to generate
% simulated voltages with a specified skip pattern. From these volatages
% fdEIT images are constructed using the mixed-minimization-algorithm.
%
% The results are saved in a mat file.

%% Specify the reconstruction and model parameters

close all
% Specify the skip pattern to be used for voltage generation and
% reconstruction
nSkip = 4;

% Specify the simulation frequencies
freq_1 = 100e3; 
freq_2 = 700e3;

% Specify thet hyper-parameters and mixed-minimization parameters to be used
%hp_param_all = exp((1:2:40))/(10^7); hp_param = hp_param_all(1:12);
hp_param=[1e-2, 5e-2, 1e-1, 5e-1:1e-1:1, 1.5, 2, 2.5, 3.5, 4, 5, 10, 20, 50];
% hp_param=[5e-1];
%hp_param=[1e-2:1e-2:9e-2, 1e-1:1e-1:9e-1, 1:1:10];
%beta=[0,0.25,0.5,0.75,1,1.25,2.5,5,7.5,10];
%beta=[0.25,0.5,0.75];
beta=[0.5];

% Specify the initial conducvitity distribution
init_iter = 0.1+1i*0.01; % initial conductivity


%% Load the pathology models 

% Inititalize the cell array to store the models and their respective names
cases = cell(1,10);

[cases{1}.fmdl,cases{1}.img_sim1,cases{1}.img_sim2,cases{1}.img_diff] = create_model_pat(freq_1, freq_2, 'healthy', 32, 0.03, nSkip);
plot_fdEITmodel(cases{1}.img_diff, 'model_healthy', 0);
cases{1}.title = 'healthy';

[cases{2}.fmdl,cases{2}.img_sim1,cases{2}.img_sim2,cases{2}.img_diff] = create_model_pat(freq_1, freq_2, 'at_right', 32, 0.03, nSkip);
plot_fdEITmodel(cases{2}.img_diff, 'model_at_right', 0);
cases{2}.title = 'at_right';

[cases{3}.fmdl,cases{3}.img_sim1,cases{3}.img_sim2,cases{3}.img_diff] = create_model_pat(freq_1, freq_2, 'at_left', 32, 0.03, nSkip);
plot_fdEITmodel(cases{3}.img_diff, 'model_at_left', 0);
cases{3}.title = 'at_left';

[cases{4}.fmdl,cases{4}.img_sim1,cases{4}.img_sim2,cases{4}.img_diff] = create_model_pat(freq_1, freq_2, 'ed_right', 32, 0.03, nSkip);
plot_fdEITmodel(cases{4}.img_diff, 'model_ed_right', 0);
cases{4}.title = 'ed_right';

[cases{5}.fmdl,cases{5}.img_sim1,cases{5}.img_sim2,cases{5}.img_diff] = create_model_pat(freq_1, freq_2, 'ed_left', 32, 0.03, nSkip);
plot_fdEITmodel(cases{5}.img_diff, 'model_ed_left', 0);
cases{5}.title = 'ed_left';

[cases{6}.fmdl,cases{6}.img_sim1,cases{6}.img_sim2,cases{6}.img_diff] = create_model_pat(freq_1, freq_2, 'peri_eff', 32, 0.03, nSkip);
plot_fdEITmodel(cases{6}.img_diff, 'model_peri_eff', 0);
cases{6}.title = 'peri_eff';

[cases{7}.fmdl,cases{7}.img_sim1,cases{7}.img_sim2,cases{7}.img_diff] = create_model_pat(freq_1, freq_2, 'ed_low_right', 32, 0.03, nSkip);
plot_fdEITmodel(cases{7}.img_diff, 'model_ed_low_right', 0);
cases{7}.title = 'ed_low_right';

[cases{8}.fmdl,cases{8}.img_sim1,cases{8}.img_sim2,cases{8}.img_diff] = create_model_pat(freq_1, freq_2, 'ed_low_left', 32, 0.03, nSkip);
plot_fdEITmodel(cases{8}.img_diff, 'model_ed_low_left', 0);
cases{8}.title = 'ed_low_left';

[cases{9}.fmdl,cases{9}.img_sim1,cases{9}.img_sim2,cases{9}.img_diff] = create_model_pat(freq_1, freq_2, 'pleu_eff_right', 32, 0.03, nSkip);
plot_fdEITmodel(cases{9}.img_diff, 'model_pleu_eff_right', 0);
cases{9}.title = 'pleu_eff_right';

[cases{10}.fmdl,cases{10}.img_sim1,cases{10}.img_sim2,cases{10}.img_diff] = create_model_pat(freq_1, freq_2, 'pleu_eff_left', 32, 0.03, nSkip);
plot_fdEITmodel(cases{10}.img_diff, 'model_pleu_eff_left', 0);
cases{10}.title = 'pleu_eff_left';


%% Generate voltage measurement

%%% Calculate the measured voltages via solving forward problem
for i = 1:length(cases)
cases{i}.v1 = fwd_solve(cases{i}.img_sim1); % The voltages at freq_1 are simulated
cases{i}.v2 = fwd_solve(cases{i}.img_sim2); % The voltages at freq_2 are simulated
end

fmdl_recon = create_thorax_fem_simon(32, nSkip, 0.05, 0, 0, 0);
figure
show_fem(mk_image(fmdl_recon,1),[0,1]);

%% Reconstruct images from the simulated voltages

for k = 8:length(cases)
    for jj = 1:length(beta)
        for ii = 1:length(hp_param)
            imdl    = my_weighted_invprob_properties(fmdl_recon, init_iter, hp_param(ii),...
                true,beta(jj));
            img_solve(ii)= inv_solve(imdl, cases{k}.v1, cases{k}.v2);
            
            disp('o-----------------------------------------------------------------------------------o')
            disp('Finished reconstructin for hyperparameter: ')
            disp(strcat(num2str(ii),':'));
            disp(num2str(hp_param(ii)));
            disp('Last reconstruction used beta number: ')
            disp(strcat(num2str(jj),':'));
            disp(num2str(beta(jj)));
            disp('Current Case Number: ')
            disp(num2str(k));
            disp('o-----------------------------------------------------------------------------------o')
            
        end
        save(strcat('recon_', cases{k}.title, '_B', num2str(jj), '_allHPs_sk44.mat'),'img_solve','beta','hp_param');
    end
    
    plot_fdEITmodel(img_solve(4), 'Healthy_Reconstruction', false);
end