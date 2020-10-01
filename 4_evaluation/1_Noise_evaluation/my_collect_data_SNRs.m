function [rms_s, rms_n, res_dB_n] = my_collect_data_SNRs

%% load data 
temp = load('RMS_s_hp851_1.mat'); rms_s_all(1,:)  = temp.RMS_s;
temp = load('RMS_n_hp851_1.mat'); rms_n_all(1,:) = temp.RMS_n;
temp = load('res_dB_n_hp851_1.mat'); res_dB_n_all(1,:) = temp.res_dB_n;

temp = load('RMS_s_hp851_2.mat'); rms_s_all(2,:)  = temp.RMS_s;
temp = load('RMS_n_hp851_2.mat'); rms_n_all(2,:) = temp.RMS_n;
temp = load('res_dB_n_hp851_2.mat'); res_dB_n_all(2,:) = temp.res_dB_n;

temp = load('RMS_s_hp851_3.mat'); rms_s_all(3,:)  = temp.RMS_s;
temp = load('RMS_n_hp851_3.mat'); rms_n_all(3,:) = temp.RMS_n;
temp = load('res_dB_n_hp851_3.mat'); res_dB_n_all(3,:) = temp.res_dB_n;

temp = load('RMS_s_hp851_4.mat'); rms_s_all(4,:)  = temp.RMS_s;
temp = load('RMS_n_hp851_4.mat'); rms_n_all(4,:) = temp.RMS_n;
temp = load('res_dB_n_hp851_4.mat'); res_dB_n_all(4,:) = temp.res_dB_n;

temp = load('RMS_s_hp851_5.mat'); rms_s_all(5,:)  = temp.RMS_s;
temp = load('RMS_n_hp851_5.mat'); rms_n_all(5,:) = temp.RMS_n;
temp = load('res_dB_n_hp851_5.mat'); res_dB_n_all(5,:) = temp.res_dB_n;


%%% Average all 5 random measurements 
rms_s=mean(rms_s_all); rms_n=mean(rms_n_all); res_dB_n=mean(res_dB_n_all); 
