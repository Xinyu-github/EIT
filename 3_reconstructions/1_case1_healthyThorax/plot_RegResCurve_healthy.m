%% Select the reconstructed parametersweep and plot the reg/res curve
close all
selection = 3;
displayAllReconstructions = 0;
displayAllHP = 0;
% name={'c1_B', 'c1_onlyGN_B', '_c1_scaling_B', '_c1_scaling2_increasedResolution_B'};
name={'c1_B', 'c1_onlyGN_B', '_c1_scaling_B', '_c1_scaling2_increasedResolution_B'};
betaNum=[8, 1, 6, 1];
% for scaling beta we can calculate opt HP based on singularity in curve
% selection = 3,4

for jj = 1:betaNum(selection)
    clear img_solve;
    % load data
    load((strcat('results/recon_',name{selection},num2str(jj),'_allHP_skip0.mat')));% here is data for '_skip0'
    [residual,regularizer] = get_norm(img_solve);
     
    %plot Reg/Res Curve
    figure(1)
    loglog(residual,regularizer);
    hold on;
    % loglog(Res,Reg);
    xlabel('complete residual');
    ylabel('regularizer');
    if displayAllHP
        for j = 1:length(img_solve)%[locPeakinv(1), locPeak, locPeakinv(2)]
        % text(residual,regularizer,num2str(labels));
            txt = strcat('\lambda = ',num2str(img_solve(j).hp));   
            text(residual(j),regularizer(j),txt,'HorizontalAlignment','right');
        end
    end
    if displayAllReconstructions
        plot_reconstructions(img_solve, 1:length(img_solve), 'Healthy Thorax');
    elseif selection >= 3
        %find interesting hp
        [pksVal, locPeak] = findpeaks(regularizer);
        [pksValinv, locPeakinv] = findpeaks(1./regularizer);
            %locPeak is peak of singularity
            %locPeakinv(1) is minimum of L-curve
            %locPeakinv(2) is first minimum after peak
        plot_reconstructions(img_solve, [locPeakinv(1), locPeak, locPeakinv(2)], 'Thorax healthy');
        
        for j = [locPeakinv(1), locPeak, locPeakinv(2)]
        % text(residual,regularizer,num2str(labels));
            figure(1)
            txt = strcat('\lambda = ',num2str(img_solve(j).hp));   
            text(residual(j),regularizer(j),txt,'HorizontalAlignment','right');
        end
    end
    if selection == 1
        betaLegend = {'beta = 0', 'beta = 0.1', 'beta = 0.25', 'beta = 0.5', ...
            'beta = 1', 'beta = 2.5', 'beta = 5', 'beta = 10'};
    elseif selection == 2
        betaLegend = {'beta = 0'};
    elseif selection == 3
        betaLegend = {'beta = 0.5*hp', 'beta = 1*hp', 'beta = 2*hp', ...
            'beta = 4*hp', 'beta = 6*hp', 'beta = 8*hp'};
    elseif selection == 4
        betaLegend = {'beta = 2*hp'};
    end
end  
figure(1);
legend(betaLegend);
title('L-curve for various beta scaling factors');