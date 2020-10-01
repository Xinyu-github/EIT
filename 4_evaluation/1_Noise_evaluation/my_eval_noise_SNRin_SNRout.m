function [res_SNR_dB_t,frac_s_n_t] = my_eval_noise_SNRin_SNRout(input_n_t, rms_s, rms_n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% input: input_noise(SNR-dB), fitobj1
%%% output: resulting_noise(SNR-dB) 
% input_n_t = 50; 
%% noise in-out relationship %     res_snr_dB = 20.*log(rms_signal/rms_noise); 
input_n = [0:5:100]; %dB

frac_s_n =  rms_s./rms_n; %(rms_signal/rms_noise)
res_SNR_dB = 20.*log(frac_s_n); 

% fitobj1 = fit(input_n',frac_s_n','smoothingspline');
% figure(10), 
% plot(fitobj1,input_n,frac_s_n), hold on,
% % plot(input_n,frac_s_n,'-k','LineWidth',2), 
% 
% xlabel('x-input noise'),ylabel('y-rmsS/rmsN')
% title('Data fitting btw input SNR, reconstructed rmsS/rmsN');

fitobj2 = fit(input_n',res_SNR_dB','smoothingspline');
figure(11),
set(gcf, 'Position', [100, 100, 1440, 700],'Color','white');
reconSNR = plot(fitobj2,input_n,res_SNR_dB), hold on,
xlabel('SNR of input signal [dB]'),ylabel('SNR of reconstruction [dB]')
%title('Data fitting btw input SNR, reconstructed SNR');

% equivalent line (input-output noise)
linetest = [0:5:100];
equivalent = plot(linetest,linetest,'Linewidth',2), hold on; %y=x

% find intersect point
[x,y]=intersections(input_n,res_SNR_dB,linetest,linetest);
disp(strcat('The intersect point btw input, reconstructed noise is ',string(x(end)),'SNR of input noise' ))
intersection = plot(x(end),y(end),'og','Linewidth',2);

legend([equivalent intersection],'equivalent line','intersection point')
grid on;
%% %%%%%%%%%%%%%%%% test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


frac_s_n_t = fitobj1(input_n_t); %input 
res_SNR_dB_t = 20*log(frac_s_n_t); 
disp(strcat('signal with SNR',string(input_n_t),...
    'dB >> reconstructed RMS_signal/RMS_noise ',string(frac_s_n_t),...
    ' >> reconstructed SNR ',string(res_SNR_dB_t),'dB'))
