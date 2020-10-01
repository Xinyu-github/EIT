function [calib_f1, calib_f2, calib_f1_pol2cart, calib_f2_pol2cart] = swissTOM_calibrate(high_ohm_mean , low_ohm_mean, data_mean_f1, data_mean_f2, show_voltages, measurement_name)
% This function calibrates swisstom measurement data using measurements
% that were taken using the calibration network that was in use at MEDIT in
% February and March of 2020. Two different calibrations are performed, one
% where only the complex numbers themselves are scaled using a calculated M
% and O and another, BETTER calibration where the measurement data is
% calibrated by magnitude and phase. This function is NOT part of the 
% EIDORS package and was developed at the MEDIT institute.
%
% Syntax: [calib_f1, calib_f2, calib_f1_pol2cart, calib_f2_pol2cart] = swissTOM_calibrate(high_ohm_mean , low_ohm_mean, data_mean_f1, data_mean_f2, show_voltages, measurement_name)
%
% Input Parameters:
% - high_ohm_mean:      Measurement of the calibration network at the
%                       higher impedance setting. This must be a matrix 
%                       of size Nx2 where N is the number of meausrements
%                       and 2 is the number of measurement frequencies. If
%                       multiple frames were measured it is advised to take
%                       a mean over the frames before giving the data to
%                       this function.
% - low_ohm_mean:       Measurement of the calibration network at the
%                       lower impedance setting. This must be a matrix 
%                       of size Nx2 where N is the number of meausrements
%                       and 2 is the number of measurement frequencies. If
%                       multiple frames were measured it is advised to take
%                       a mean over the frames before giving the data to
%                       this function.
% - data_mean_f1:       Input measurement vector at lower measurement
%                       frequency. If multiple frames were measured it is 
%                       advised to take a mean over the frames before 
%                       giving the data to this function.
% - data_mean_f2:       Input measurement vector at higher measurement
%                       frequency. If multiple frames were measured it is 
%                       advised to take a mean over the frames before 
%                       giving the data to this function.
% - show_voltages:      Parameter that specifies if the voltages before and
%                       after calibration shall be plotted and displayed.
% - measurement_name:   Name of the measurement that is currently being
%                       calibrated. This string is used as a title for the
%                       voltage plots mentioned in the parameter
%                       description above.
%
% Output Parameters:
% - calib_f1:           Calibrated data at lower measurement frequency. The
%                       calibration was performed only by scaling the
%                       complex raw numbers and applying an offset to them.
% - calib_f2:           Calibrated data at higher measurement frequency. The
%                       calibration was performed only by scaling the
%                       complex raw numbers and applying an offset to them
% - calib_f1_pol2cart:  Calibrated data at lower frequency that was
%                       calibrated using the advanced calibration by
%                       magnitude and phase.
% - calib_f2_pol2cart:  Calibrated data at higher frequency that was
%                       calibrated using the advanced calibration by
%                       magnitude and phase.

% Calculate gradient and offset from the measurements taken with the
% calibration network
M = ((0.1-8)./(abs(low_ohm_mean)-abs(high_ohm_mean)));
O = 8 - M.*abs(high_ohm_mean);

% Display the voltages before calibration. The magnitude, real 
% and imag of the measured voltages is displayed.
if show_voltages
    
    xaxis = [1:1:928];
    
    figure;
    subplot(3,1,1)
    plot(abs(data_mean_f1))
    axis([1, 928, min(abs([data_mean_f1; data_mean_f2])), max(abs([data_mean_f1; data_mean_f2]))])
    title(strcat('Abs of voltages before calibration (f_low), ', measurement_name), 'Interpreter', 'none')
    subplot(3,1,2)
    plot(abs(data_mean_f2))
    axis([1, 928, min(abs([data_mean_f1; data_mean_f2])), max(abs([data_mean_f1; data_mean_f2]))])
    title(strcat('Abs of voltages before calibration (f_high), ', measurement_name), 'Interpreter', 'none')
    subplot(3,1,3)
    plot(abs(data_mean_f2))
    title(strcat('Abs of voltages difference before calibration, ', measurement_name), 'Interpreter', 'none')
    xlim([1, 928])
    
    figure;
    subplot(3,2,1)
    plot(real(data_mean_f1))
    title(strcat('Real of voltages before calibration (f_low), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(real([data_mean_f1; data_mean_f2])), max(real([data_mean_f1; data_mean_f2]))])
    subplot(3,2,3)
    plot(real(data_mean_f2))
    title(strcat('Real of voltages before calibration (f_high), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(real([data_mean_f1; data_mean_f2])), max(real([data_mean_f1; data_mean_f2]))])
    subplot(3,2,5)
    plot(real(data_mean_f2))
    title(strcat('Real of voltages difference before calibration, ', measurement_name), 'Interpreter', 'none')
    xlim([1, 928])
    
    subplot(3,2,2)
    plot(imag(data_mean_f1))
    title(strcat('Im of voltages before calibration (f_low), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(imag([data_mean_f1; data_mean_f2])), max(imag([data_mean_f1; data_mean_f2]))])
    subplot(3,2,4)
    plot(imag(data_mean_f2))
    title(strcat('Im of voltages before calibration (f_high), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(imag([data_mean_f1; data_mean_f2])), max(imag([data_mean_f1; data_mean_f2]))])
    subplot(3,2,6)
    plot(imag(data_mean_f2))
    title(strcat('Im of voltages difference before calibration, ', measurement_name), 'Interpreter', 'none')
    xlim([1, 928])
    
end

calib_0_1Ohm_mean = low_ohm_mean;

% Calibration by magnitude and phase
calib_abs_f1 = (abs(data_mean_f1).*M(:, 1) + O(:, 1));
calib_angle_f1 = angle(data_mean_f1) - angle(calib_0_1Ohm_mean(:, 1));
calib_abs_f2 = (abs(data_mean_f2).*M(:, 2) + O(:, 2));
calib_angle_f2 = angle(data_mean_f2) - angle(calib_0_1Ohm_mean(:, 2));

[calib_f1_x, calib_f1_y] = pol2cart(calib_angle_f1, calib_abs_f1);
calib_f1_pol2cart = complex(calib_f1_x, calib_f1_y);
[calib_f2_x, calib_f2_y] = pol2cart(calib_angle_f2, calib_abs_f2);
calib_f2_pol2cart = complex(calib_f2_x, calib_f2_y);

% Simple, direct scaling of cartesian complex numbers
calib_f1 = (data_mean_f1.*M(:, 1) + O(:, 1))./32;
calib_f2 = (data_mean_f2.*M(:, 2) + O(:, 2))./32;

% Display the voltages after calibration. The magnitude, real 
% and imag of the measured voltages is displayed.
if show_voltages
    figure;
    subplot(3,1,1)
    plot(abs(calib_f1_pol2cart))
    title(strcat('Abs of voltages after calibration (f_low), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(abs([calib_f1_pol2cart; calib_f2_pol2cart])), max(abs([calib_f1_pol2cart; calib_f2_pol2cart]))])
    subplot(3,1,2)
    plot(abs(calib_f2_pol2cart))
    title(strcat('Abs of voltages after calibration (f_high), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(abs([calib_f1_pol2cart; calib_f2_pol2cart])), max(abs([calib_f1_pol2cart; calib_f2_pol2cart]))])
    subplot(3,1,3)
    plot(abs(calib_f2_pol2cart-calib_f1_pol2cart))
    title(strcat('Abs of voltages difference after calibration, ', measurement_name), 'Interpreter', 'none')
    xlim([1, 928])
    
    figure;
    subplot(3,2,1)
    plot(real(calib_f1_pol2cart))
    title(strcat('Real of voltages after calibration (f_low), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(real([calib_f1_pol2cart; calib_f2_pol2cart])), max(real([calib_f1_pol2cart; calib_f2_pol2cart]))])
    subplot(3,2,3)
    plot(real(calib_f2_pol2cart))
    title(strcat('Real of voltages after calibration (f_high), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(real([calib_f1_pol2cart; calib_f2_pol2cart])), max(real([calib_f1_pol2cart; calib_f2_pol2cart]))])
    subplot(3,2,5)
    plot(real(calib_f2_pol2cart-calib_f1_pol2cart))
    title(strcat('Real of voltages difference after calibration, ', measurement_name), 'Interpreter', 'none')
    xlim([1, 928])
    
    subplot(3,2,2)
    plot(imag(calib_f1_pol2cart))
    title(strcat('Im of voltages after calibration (f_low), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(real([calib_f1_pol2cart; calib_f2_pol2cart])), max(real([calib_f1_pol2cart; calib_f2_pol2cart]))])
    subplot(3,2,4)
    plot(imag(calib_f2_pol2cart))
    title(strcat('Im of voltages after calibration (f_high), ', measurement_name), 'Interpreter', 'none')
    axis([1, 928, min(real([calib_f1_pol2cart; calib_f2_pol2cart])), max(real([calib_f1_pol2cart; calib_f2_pol2cart]))])
    subplot(3,2,6)
    plot(imag(calib_f2_pol2cart-calib_f1_pol2cart))
    title(strcat('Img of voltages difference after calibration, ', measurement_name), 'Interpreter', 'none')
    xlim([1, 928])
    
end

end

