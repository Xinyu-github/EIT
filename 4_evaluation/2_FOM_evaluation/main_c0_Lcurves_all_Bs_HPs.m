% clc,clear;
%hp_param_all = exp((0.4:8:100)/2)/(5*10^19); hp_param = hp_param_all(3:11);
%B=[0.01, 0.1, 0.25, 1];
hp_param_all = exp((0.4:8:100)/2)/(5*10^19); hp_param = hp_param_all(4:7);
B=[0.01, 0.1, 0.25, 1];

x_term = []; y_term = []; hp = []; ls=[]; mixed=[];

%% Load data from the saved files (Reconstruction solutions) and Generate L-curves
for j=1:length(B)
    reconName = strcat('recon_FOM_B',num2str(B(j)),'cyl2cyl');
    reconName(regexp(reconName,'[.]'))=[];
    eq = load(reconName); eq = eq.eq_struct;
    [x_term(j,:), y_term(j,:), hp(j,:), ls(j,:), mixed(j,:)] = my_Lcurve_c2_1(eq,B(j));
    %figure(100+j);plot(x_term(j,:), y_term(j,:),'-','linewidth',2),hold on,
end

%% Select acceptable B whose absolute-mixed minimisations do not dominate the reconstruction
%%% in case: user want absolute term be less than 2*(residual term) or else
ls_base = 2.*max(ls(:));% magnitude of least-square/residual term
figure(102), plot(hp(1,:), ones(1,size(mixed,2)).*ls_base,    'ok', 'linewidth', 2), hold on,
for i=1:size(mixed,1) %plot all magnitudes
     plot(hp(i,:), mixed(i,:),     'o', 'linewidth', 2), hold on,
end

mixed_base = []; % magnitude of absolute-mixed minimising terms
for b=1:size(mixed,1) %for each B
    mixed_base(b,1) = mean(mixed(b,:)); % magnitude of absolute-mixed minimising terms
end
selected_Bs_indx = find(mixed_base(:) < ls_base); % Acceptable Bs


%% Calculate Curvatures for all selected Bs
HPindx_maxK = []; % collects max curvature of each L-curve (at corner)
for i=1:length(selected_Bs_indx)%i is number of beta
dx  = gradient(x_term(i,:));
ddx = gradient(dx);
dy  = gradient(y_term(i,:));
ddy = gradient(dy);
num   = dx .* ddy - ddx .* dy;
denom = dx .* dx + dy .* dy;
denom = sqrt(denom);
denom = denom.^3; %* denom * denom;
curvature = num ./ denom;
curvature(denom < 0) = NaN;
curvature = max(curvature)-curvature;

%%% Sort indices to plot
Kcurv(i,:) = curvature(end:-1:1);
HP(i,:) = hp(i,end:-1:1);

figure(101), plot(HP(i,:), Kcurv(i,:), 'Linewidth', 2), hold on,
end

figure(100), xlabel('Residual term with absolute-mixed minimisation'), 
ylabel('Regularisation term'), axis([0 1 0 1]),
legend('\beta=0.25','\beta=0.5','\beta=0.75','\beta=1.00','\beta=1.25','\beta=2.50','\beta=5.00'),
title('4.1) Normalised L-curves for all variant \lambdas or HPs'), 


figure(102), legend('mixed term \beta=0.25','mixed term \beta=0.5','mixed term \beta=0.75',...
'mixed term \beta=1.00','mixed term \beta=1.25','mixed term \beta=2.50','mixed term \beta=5.00');
title('4.2) Comparison of residual term and mixed terms for different \beta s'), 
ylabel('Magnitude of terms'), xlabel('hyper-parameter'), grid on;


figure(101), legend('\beta=0.25','\beta=0.5','\beta=0.75','\beta=1.00','\beta=1.25');
title('4.3) Curvatures for all selected \betas'), ylabel('Calculated Curvature'), xlabel('hyper-parameter'), grid on;

