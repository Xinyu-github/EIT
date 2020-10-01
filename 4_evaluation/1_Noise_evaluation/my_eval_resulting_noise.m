function [res_SNR_in_dB,rms_signal,rms_noise,res_img_noise] = my_eval_resulting_noise(res_img_Lower_in_dB,res_img_Higher_in_dB)

%% resulting noise subtracted by resulting images of no noise and SNR ....
img1 = res_img_Lower_in_dB; %img_solve_w2_0;
img2 = res_img_Higher_in_dB;%img_solve_w2_40;
img3 = img1;

res_noise = img2.elem_data(:) - img1.elem_data(:);
img3.elem_data = res_noise;

sum_noise = sum(res_noise(:));
abs_noise = abs(sum_noise);
max_noise = max(res_noise);

%% calculate SNR SNR(dB) = 10log(Psignal/Pnoise) = 10log(RMSsignal^2 / RMSnoise^2)
rms_noise = rms(res_noise);
rms_signal = rms(img1.elem_data(:));

res_SNR_in_dB = 20*log(rms_signal) - 20*log(rms_noise);
res_img_noise = img3;