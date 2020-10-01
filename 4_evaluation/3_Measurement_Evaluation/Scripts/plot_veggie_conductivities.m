% This script can be used to display the measured vegetable conductivities
% that are stored in the folder '2020_03_05_veggie_conductivities'. All
% relevant parameters can be set below.

%% Set parameters
measurement_style = 'short'; % Specify the measurement type that shall be display ('short' or 'long')
no_of_measurements = 60; % Specify if the measurement with 30 or 60 datapoints shall be displayed
onlye4toe6 = true; % Specify if ONLY the measurements taken between 10kHz and 1MHz shall be displayed

%% Load specified data

measurements = load('.\7_Measurement_Evaluation\2020_03_05_veggie_conductivities\conductivity_measurements_vegetables.mat');

% Load the right data according to the specifications above
if strcmp(measurement_style, 'long') && (no_of_measurements == 60)
    beep
    disp('There are no measurements for this combination of parameters.');
else
    % Store all measurements in a struct (Measurements are impedances)
    imp_pot_e3_e5 = measurements.(strcat('BIS_pot_', measurement_style, '_e5_', num2str(no_of_measurements)));
    imp_pot_e3_e6 = measurements.(strcat('BIS_pot_', measurement_style, '_e6_', num2str(no_of_measurements)));
    imp_pot_e4_e5 = measurements.(strcat('BIS_pot_', measurement_style, '_e4_e5_', num2str(no_of_measurements)));
    imp_pot_e4_e6 = measurements.(strcat('BIS_pot_', measurement_style, '_e4_e6_', num2str(no_of_measurements)));
    
    imp_pump_e3_e5 = measurements.(strcat('BIS_pump_', measurement_style, '_e5_', num2str(no_of_measurements)));
    imp_pump_e3_e6 = measurements.(strcat('BIS_pump_', measurement_style, '_e6_', num2str(no_of_measurements)));
    imp_pump_e4_e5 = measurements.(strcat('BIS_pump_', measurement_style, '_e4_e5_', num2str(no_of_measurements)));
    imp_pump_e4_e6 = measurements.(strcat('BIS_pump_', measurement_style, '_e4_e6_', num2str(no_of_measurements)));
end

% Invert the measurements since they were measured in an inverted way
imp_pot_e4_e5 = imp_pot_e4_e5(end:-1:1);
imp_pot_e4_e6 = imp_pot_e4_e6(end:-1:1);
imp_pump_e4_e5 = imp_pump_e4_e5(end:-1:1);
imp_pump_e4_e6 = imp_pump_e4_e6(end:-1:1);

% Create the x-axes
x_e3_e5 = 1e3:(1e5-1e3)/(no_of_measurements-1):1e5;
x_e3_e6 = 1e3:(1e6-1e3)/(no_of_measurements-1):1e6;
x_e4_e5 = 1e4:(1e5-1e4)/(no_of_measurements-1):1e5;
x_e4_e6 = 1e4:(1e6-1e4)/(no_of_measurements-1):1e6;

% Convert to conductivities
cond_pot_e3_e5 = 1./imp_pot_e3_e5;
cond_pot_e3_e6 = 1./imp_pot_e3_e6;
cond_pot_e4_e5 = 1./imp_pot_e4_e5;
cond_pot_e4_e6 = 1./imp_pot_e4_e6;

cond_pump_e3_e5 = 1./imp_pump_e3_e5;
cond_pump_e3_e6 = 1./imp_pump_e3_e6;
cond_pump_e4_e5 = 1./imp_pump_e4_e5;
cond_pump_e4_e6 = 1./imp_pump_e4_e6;

%% Plot data

if onlye4toe6
    
    figure('Name', strcat('Measured Potato and Pumpkin Conductivities [', measurement_style, ' measurement, ', num2str(no_of_measurements), ' data points]'));
    
    subplot(2,2,1)
    semilogx(x_e4_e6, abs(cond_pot_e4_e6)*(0.028/(0.037*0.034)), 'LineWidth', 2);
    title('Kartfoffel');
    xlabel('Frequenz [Hz]');
    ylabel('Betrag [S/m]');
    set(gca,'FontName','Times')
    
    subplot(2,2,2)
    semilogx(x_e4_e6, rad2deg(angle(cond_pot_e4_e6)), 'LineWidth', 2);
    title('Kartfoffel');
    xlabel('Frequenz [Hz]');
    ylabel('Phase [°]');
    set(gca,'FontName','Times')
    
    subplot(2,2,3)
    semilogx(x_e4_e6, abs(cond_pump_e4_e6)*(0.028/(0.037*0.034)), 'LineWidth', 2);
    ylim([0, 0.25]);
    title('Kürbis');
    xlabel('Frequenz [Hz]');
    ylabel('Betrag [S/m]');
    set(gca,'FontName','Times')
    
    subplot(2,2,4)
    semilogx(x_e4_e6, rad2deg(angle(cond_pump_e4_e6)), 'LineWidth', 2);
    title('Kürbis');
    xlabel('Frequenz [Hz]');
    ylabel('Phase [°]');
    set(gca,'FontName','Times')
    
    % Calculate the differential conductivities between 77,12kHz and
    % 194,58khz
    complex_diff_pot = (cond_pot_e4_e6(12) - cond_pot_e4_e6(5)).*(0.028/(0.037*0.034));
    complex_diff_pump = (cond_pump_e4_e6(12) - cond_pump_e4_e6(5)).*(0.028/(0.037*0.034));
    
    % Display the differentials in numbers
    disp('Magnitude, Phase and Complex of the differential conductivity of the potato (77119Hz, 194580Hz):');
    disp(abs(complex_diff_pot));
    disp(rad2deg(angle(complex_diff_pot)));
    disp(complex_diff_pot);
    disp('Magnitude, Phase and Complex of the differential conductivity of the pumpkin (77119Hz, 194580Hz):');
    disp(abs(complex_diff_pump));
    disp(rad2deg(angle(complex_diff_pump)));
    disp(complex_diff_pump);
    
    % Plot the differentials in the imaginary plane
    figure
    hold all
    grid on
    p1 = plot(cond_pot_e4_e6.*(0.028/(0.037*0.034)), 'LineWidth', 2);
    p2 = plot(cond_pump_e4_e6.*(0.028/(0.037*0.034)), 'LineWidth', 2);
    
    plot(cond_pump_e4_e6(5)*(0.028/(0.037*0.034)), 'o', 'color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'k')
    plot(cond_pump_e4_e6(12)*(0.028/(0.037*0.034)), 'o', 'color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'k')
    plot(cond_pot_e4_e6(12)*(0.028/(0.037*0.034)), 'o', 'color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'k')
    plot(cond_pot_e4_e6(5)*(0.028/(0.037*0.034)), 'o', 'color', 'k', 'MarkerSize', 5, 'MarkerFaceColor', 'k')
    
    title('Complex conductivities from 10kHz to 1MHz (Markers: 77119Hz, 194580Hz)');
    xlabel('Real\{\gamma(\omega)\} [S/m]');
    ylabel('Imag\{\gamma(\omega)\} [S/m]');
    legend([p1, p2], {'Kartoffel', 'Kürbis'});
    set(gca,'FontName','Times')
    
    hold off
    
else
    % Plot all measured conductivities
    
    % Plot potato data
    figure('Name', strcat('Measured Potato Conductivities [', measurement_style, ' measurement, ', num2str(no_of_measurements), ' data points]'));
    
    subplot(4,2,1);
    semilogx(x_e3_e5, abs(cond_pot_e3_e5)*(0.028/(0.037*0.034)));
    title('Magnitude Potato 1kHz bis 100kHz');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [S/m]');
    
    subplot(4,2,2);
    semilogx(x_e3_e5, rad2deg(angle(cond_pot_e3_e5)));
    title('Phase Potato 1kHz bis 100kHz');
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    
    subplot(4,2,3);
    semilogx(x_e3_e6, abs(cond_pot_e3_e6)*(0.028/(0.037*0.034)));
    title('Magnitude Potato 1kHz bis 1MHz');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [S/m]');
    
    subplot(4,2,4);
    semilogx(x_e3_e6, rad2deg(angle(cond_pot_e3_e6)));
    title('Phase Potato 1kHz bis 1MHz');
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    
    subplot(4,2,5);
    semilogx(x_e4_e5, abs(cond_pot_e4_e5)*(0.028/(0.037*0.034)));
    title('Magnitude Potato 10kHz bis 100kHz');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [S/m]');
    
    subplot(4,2,6);
    semilogx(x_e4_e5, rad2deg(angle(cond_pot_e4_e5)));
    title('Phase Potato 10kHz bis 100kHz');
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    
    subplot(4,2,7);
    semilogx(x_e4_e6, abs(cond_pot_e4_e6)*(0.028/(0.037*0.034)));
    title('Magnitude Potato 10kHz bis 1MHz');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [S/m]');
    
    subplot(4,2,8);
    semilogx(x_e4_e6, rad2deg(angle(cond_pot_e4_e6)));
    title('Phase Potato 10kHz bis 1MHz');
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    
    % Plot pumpkin data
    figure('Name', strcat('Measured Pumpkin Conductivities [', measurement_style, ' measurement, ', num2str(no_of_measurements), ' data points]'));
    
    subplot(4,2,1);
    semilogx(x_e3_e5, abs(cond_pump_e3_e5)*(0.028/(0.037*0.034)));
    title('Magnitude Pumkpin 1kHz bis 100kHz');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [S/m]');
    
    subplot(4,2,2);
    semilogx(x_e3_e5, rad2deg(angle(cond_pump_e3_e5)));
    title('Phase Pumpkin 1kHz bis 100kHz');
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    
    subplot(4,2,3);
    semilogx(x_e3_e6, abs(cond_pump_e3_e6)*(0.028/(0.037*0.034)));
    title('Magnitude Pumpkin 1kHz bis 1MHz');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [S/m]');
    
    subplot(4,2,4);
    semilogx(x_e3_e6, rad2deg(angle(cond_pump_e3_e6)));
    title('Phase Pumpkin 1kHz bis 1MHz');
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    
    subplot(4,2,5);
    semilogx(x_e4_e5, abs(cond_pump_e4_e5)*(0.028/(0.037*0.034)));
    title('Magnitude Pumpkin 10kHz bis 100kHz');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [S/m]');
    
    subplot(4,2,6);
    semilogx(x_e4_e5, rad2deg(angle(cond_pump_e4_e5)));
    title('Phase Pumpkin 10kHz bis 100kHz');
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    
    subplot(4,2,7);
    semilogx(x_e4_e6, abs(cond_pump_e4_e6)*(0.028/(0.037*0.034)));
    title('Magnitude Pumpkin 10kHz bis 1MHz');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [S/m]');
    
    subplot(4,2,8);
    semilogx(x_e4_e6, rad2deg(angle(cond_pump_e4_e6)));
    title('Phase Pumpkin 10kHz bis 1MHz');
    xlabel('Frequency [Hz]');
    ylabel('Phase [deg]');
    
end