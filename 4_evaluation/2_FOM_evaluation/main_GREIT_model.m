%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% GREIT Test Parameters %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Reconstruct the input model with the proposed reconstruction algorithm
figure(1), set(gcf,'name','Input model for test GREIT','numbertitle','off')
%figure(300),set(gcf,'name','Reconstruction (Absolute)','numbertitle','off')

% conductivity diff lung: abs(0.0212 + 0.0196i) = 0.0289
% conductivity diff back: abs(0.1124 + 0.0593i) = 0.1271

%% Create models with targets at different positions r
r = 0:0.1:0.8;
for k=1:length(r)
    cyl_center = [r(k),0];
    [v1, v2, v1w, v2w, alpha1, imdl, img_diff] = create_cyl_model(cyl_center);
    img_diff_show = img_diff; img_diff_show.elem_data = abs(img_diff.elem_data);
    figure(1), subplot(3,3,k), show_slices(img_diff_show);eidors_colourbar(img_diff_show);
    vLw(k,:) = v1w.meas;
    vHw(k,:) = v2w.meas;    
    alpha1w(k) = alpha1;
    vL(k,:) = v1.meas;
    vH(k,:) = v2.meas; 
end

%% Reconstruct with various HPs
imdl_fine = mk_common_model('i2C'); %generate fwd_model for reconstruction
hp_param=logspace(-4,1.5,30);
init_iter = 0.01+1i*0.01;%0.01+1i*0.01; 
img_solve = reconstruct_multipleHP(imdl_fine.fwd_model, hp_param, 2*hp_param, init_iter, vL(1,:)', vH(1,:)', 'FOM_fdEIT.mat');


%background: 0.1124 + 0.0593i
%Kulthisa Greit: 0.005+1i*0.005
%Kulthisa Thorax: 0.01+1i*0.01;
%Sebastien: 0.5+1i*0.03
%Here - 0.01+1i*0.01

%% fd-Reconstruction
figure,
init_iter = 0.01+1i*0.01;
hp_param_all = exp((0.4:8:100)/2)/(5*10^19);
hp1 = hp_param_all(9);% 6-for beta 0.01; 7-for beta 0.1
beta = 0.1;
imdl2    = my_weighted_invprob_properties(imdl, init_iter, hp1, alpha1, beta);
for k=1:length(r)
    imgS{k} = inv_solve(imdl2, vLw(k,:)', vHw(k,:)');
    subplot(3,3,k),
    imgS_show = imgS{k};
    imgS_show.elem_data = abs(imgS{k}.elem_data);
    show_slices(imgS_show);
    eidors_colourbar(imgS_show)
end

save('all_recons_GREIT_r03_beta01_cyl-h2c.mat', 'imgS');
%imgS = load('all_recons_GREIT_r03_beta01.mat'); imgS = imgS.imgS;

%% GN
%find the best hp
% imdl.fwd_model = mdl_normalize(imdl.fwd_model, 1); %1-only for difference (amplitude is more precise)
% hp_param_all = exp((0.4:8:100)/2)/(5*10^19); 
% k=1
% figure;
% for i=1:9
%     imdl.hyperparameter.value=hp_param_all(i+4)
%     imgGN = inv_solve(imdl, vLw(k,:)', vHw(k,:)');
%     subplot(3,3,i),
%     imgS_show = imgGN;
%     imgS_show.elem_data = abs(imgGN.elem_data);
%     show_slices(imgS_show);
%     eidors_colourbar(imgS_show)
% end
%9-best for GN
figure;
hp_param_all = exp((0.4:8:100)/2)/(5*10^19); 
imdl.fwd_model = mdl_normalize(imdl.fwd_model, 1); %1-(amplitude is more precise)
imdl.hyperparameter.value= hp_param_all(9);
for k=1:length(r)
    imgGN{k} = inv_solve(imdl, vLw(k,:)', vHw(k,:)');
    subplot(3,3,k),
    imgS_show = imgGN{k};
    imgS_show.elem_data = abs(imgGN{k}.elem_data);
    show_slices(imgS_show);
    eidors_colourbar(imgS_show)
end

%% GREIT
%imdl_greit = mk_common_model('f3cr',16);
fmdl_greit = ng_mk_cyl_models([2 1 0.05] ,[16,1],[0.1]);
fmdl_greit.stimulation =  mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
fmdl_greit = mdl_normalize(fmdl_greit, 1); %1- more precise amplitude, but worse AR
%imdl_greit = mdl_normalize(imdl_greit.fwd_model, 1);
% img = mk_image(fmdl_greit,1); % Homogeneous background
% r =  linspace(0,0.9,9);
% xyzr = [r; zeros(1,9); ones(1,9);
%      0.05*ones(1,9)];
% 
% [vh,vi] = simulate_movement(img, xyzr);

opt.noise_figure = 0.5; 
opt.distr = 0; % best for cylinders
i_grc = mk_GREIT_model(fmdl_greit,.2,[],opt); 
%i_grc.normalize_measurements = 1;
figure;
for k=1:length(r)
    imgGREIT{k} = inv_solve(i_grc, vLw(k,:)', vHw(k,:)');
    %imgGREIT{k} = inv_solve(i_grc, vh, vi);
    subplot(3,3,k),
    imgS_show = imgGREIT{k};
    imgS_show.elem_data = abs(imgGREIT{k}.elem_data);
    show_slices(imgS_show);
    eidors_colourbar(imgS_show)
end

%% Backprojection
figure;
i_bp = mk_common_gridmdl('backproj');
for k=1:length(r)
    imgBP{k} = inv_solve(i_bp, vLw(k,:)', vHw(k,:)');
    subplot(3,3,k),
    imgS_show = imgBP{k};
    imgS_show.elem_data = abs(imgBP{k}.elem_data);
    show_slices(imgS_show);
    eidors_colourbar(imgS_show)
end

%save('compare_algorithms_GREIT_r01_beta01_cyl-h2c.mat', 'imgS', 'imgBP', 'imgGN', 'imgGREIT');

    %load 'colormaps.mat'
    
    %colormap default
    
%calc_colours('backgnd',[.8,.8,.8]);    % gray background colour
%calc_colours('cmap_type','draeger');  % Draegerwerk colourmap

showRecon = 0;
showComparison = 0;
showFigMerit = 1;
fSize = 18;

if showFigMerit == 1
    figure(500);
    set(gcf, 'Position', [100, 100, 1440, 700],'Color','white');
end

if showComparison == 1
    figure(400);
    set(gcf, 'Position', [100, 100, 1440, 700],'Color','white');
end


for j=1:4
    if j==1
        imgTemp = imgBP;
    elseif j==2
        imgTemp = imgGN;
    elseif j==3
        imgTemp = imgGREIT;
    else
        imgTemp = imgS;
    end 
    
    if showRecon ==1
        for k=1:length(r)    
            img_show = imgTemp{k};
            img_show.elem_data = abs(imgTemp{k}.elem_data);
            figure(300+j), subplot(3,3,k), show_slices(img_show);
            eidors_colourbar(img_show);
        end
    end
    if showComparison == 1
        numPlot=[1,5,9];
        for k=1:3    
            img_show = imgTemp{numPlot(k)};
            img_show.elem_data = abs(imgTemp{numPlot(k)}.elem_data);
            if j<3 %for better axes scaling
                img_show.elem_data = img_show.elem_data*1000;
            end
            if j~=4
                img_show.calc_colours.ref_level = max(img_show.elem_data)/2;
            end
            figure(400), subplot(3,4,(k-1)*4+j), out_img = show_slices(img_show);
            eidors_colourbar(img_show);
        end 
    end

    if showFigMerit == 1
            imgr = imgTemp{1};  
        for k=2:length(r)
            imgr.elem_data(:,k) = imgTemp{k}.elem_data;
        end



        %% GREIT test parameters
        % calculate the GREIT parameters
        R = 0.1;%cyl_radius; % target radius
        Npos = 9; % number of positions
        Xpos =  r; % positions to simulated along x-axis
        Ypos = zeros(1,Npos); 
        xyr = [Xpos; Ypos; R*ones(1,Npos)];

        Zpos = ones(1,Npos);  %% for off-plane, adjust the level (*1.5)
        xyzr = [Xpos; Ypos; Zpos; R*ones(1,Npos)];

        imgr.elem_data = abs(imgr.elem_data);
        %levels =[inf,inf,1];
        %figure(2),show_slices(imgr, levels);
        %imgr.calc_colours.npoints = 128;
        %imgr.calc_slices.levels=levels;
        if j~=4
            params = eval_GREIT_fig_merit(imgr, xyr);
        else
            params = eval_GREIT_fig_merit_fdEIT(imgr, xyr);
        end

         %RNG is not applicable due to only positive absolute freq-diff values
        % p_names = {'AR','PE','RES','SD','RNG'};
        % for i=1:5; subplot(5,1,i);
        %     plot(params(i,:)); ylabel(p_names{i});
        % end
        figure(500);
        subplot(4,1,1),
        ar{j}=plot(Xpos, params(1,:)/max(params(1,:))); 
        axis([0,Xpos(end),0.1, 1.1]);
        set(gca,'FontSize',fSize),hold on,ylabel('AR'),grid on,
        subplot(4,1,2),
        pe{j}=plot(Xpos, params(2,:)); 
        axis([0,Xpos(end),-0.05, 0.15]);  
        set(gca,'FontSize',fSize),hold on,ylabel('PE'),grid on,
        subplot(4,1,3),
        res=plot(Xpos, params(3,:)); 
        axis([0,Xpos(end),0.1, 0.4]);  
        set(gca,'FontSize',fSize),hold on,ylabel('RES'),grid on,
        subplot(4,1,4),
        sd=plot(Xpos, params(4,:)); 
        axis([0,Xpos(end),-0.05, 0.25]);  
        set(gca,'FontSize',fSize),hold on,ylabel('SD'),grid on,xlabel('target position of inhomogenity away from center');

    %     subplot(5,1,5), plot(Xpos, params(5,:)); 
    %     axis([0,Xpos(end),-2.0, 2.0]);  
    %     ylabel('RNG'), xlabel('CYL target position away to center');
    end
end
if showFigMerit == 1
    figure(500);
        hL = legend([pe{1} pe{2} pe{3} pe{4}],{'Backprojection','Gauss-Newton','GREIT','Mixed-Minimization'});
        newPosition = [0.78 0.68 0.14 0.14];
        newUnits = 'normalized';
        set(hL,'Position', newPosition,'Units', newUnits,'FontSize', 18);
        %set(hL,'Location', 'southeast','FontSize', 14);
    export_fig 'result_fig_of_merit' -pdf
end

if showComparison == 1
    figure(400);
    load 'colormaps.mat'
    cm_0.abs(1,:) = [0.6, 0.6, 0.6];
    colormap(cm_0.abs)
    %export to corel -> pdf xxxx export_fig 'result_fig_of_merit_recon' -pdf
end