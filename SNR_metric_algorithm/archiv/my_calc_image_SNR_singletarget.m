function [SNRmean,imgrs, SE, debug] = my_calc_image_SNR_singletarget(imdl, hyperparameter,find_min, doPlot)
%% CALC_IMAGE_SNR: Calculates the signal-to-noise ratio (SNR) in the image 
% domain as proposed by Braun et al. in:
% F Braun et al., A Versatile Noise Performance Metric for Electrical
% Impedance Tomography Algorithms, IEEE Trans. Biomed. Eng. 2017 (submitted).
%
%   [SNRmean, SE, debug] = calc_image_SNR(imdl, hyperparameter, doPlot)
%
% Output:
%   SNRmean   - mean of all SNR values as in equ. (9) in publication below
%   SE        - std of all SNR values as in equ. (9) in publication below
%   debug     - structure holding some information used for debug purposes
%
% Input:
%   imdl            - inverse model (EIDORS struct)
%      imdl.hyperparameter.roi_scaling_factor - amount of model shrinking
%                                   to determine ROI where seed n_T targets
%                                   (DEFAULT: 0.5);
%      imdl.hyperparameter.n_targets - number of targets to seed in ROI n_T
%                                   (DEFAULT: n_T = 200);
%      imdl.hyperparameter.target_radius - relative target radius r_T
%                                   (DEFAULT: r_T = 0.05);
%      imdl.hyperparameter.xyzr_targets - vector 4 x num targets n_T
%                                   specify targets manually [x y z r]
%   hyperparameter  - desired hyperparameter value, this will overwrite the
%                     imdl.hyperparameter.value and is for compatibility
%                     purposes with the function calc_noise_figure()
%   doPlot    will enable plotting if set to true (default = false)
%
%
% See also: CALC_NOISE_FIGURE
%
% Fabian Braun, December 2016
%
% CITATION_REQUEST:
% AUTHOR: F Braun et al.
% TITLE: A Versatile Noise Performance Metric for Electrical Impedance Tomography Algorithms
% JOURNAL: IEEE Transactions on Biomedical Engineering
% YEAR: 2017
% VOL: 64
% NUM: 10
% PAGE: 2321-2330
% DOI: 10.1109/TBME.2017.2659540
%

% (C) 2016 Fabian Braun. License: GPL version 2 or version 3
% $Id: calc_image_SNR.m 5646 2017-09-22 06:16:33Z fab-b $

citeme(mfilename);

%% Configuration
% model shrunken by this factor delimits the ROI
ROI_SCALING_FACTOR = 0.5; 
% approximate number of targets to uniformly distribute in the ROI
N_DESIRED_TARGETS = 200;
% relative radius of the circular(2D)/spherical(3D) target (this is relative to the outer model radius) 
NORM_TGT_RADIUS = 0.05; % original
% TEZ region threshold 
TEZ_REGION_THRESH = 0.25;  % as ratio of maximum value : QUATER AMPLITUDE THRESHOLD


%% Parse and default inputs
if ~exist('doPlot', 'var') || isempty(doPlot)
    doPlot = false;
end
if isfield(imdl.hyperparameter, 'roi_scaling_factor') && ~isempty(imdl.hyperparameter.roi_scaling_factor)
    ROI_SCALING_FACTOR = imdl.hyperparameter.roi_scaling_factor;
end
assert(ROI_SCALING_FACTOR < 1, 'ROI must be scaled smaller than effective model');
if isfield(imdl.hyperparameter, 'n_targets') && ~isempty(imdl.hyperparameter.n_targets)
    N_DESIRED_TARGETS = imdl.hyperparameter.n_targets;
end
if isfield(imdl.hyperparameter, 'target_radius') && ~isempty(imdl.hyperparameter.target_radius)
    NORM_TGT_RADIUS = imdl.hyperparameter.target_radius;
end
if isfield(imdl.hyperparameter, 'xyzr_targets') && ~isempty(imdl.hyperparameter.xyzr_targets)
    xyzr_targets = imdl.hyperparameter.xyzr_targets;
else
    xyzr_targets = [];
end


%% input parsing
if nargin>=2 && numel(hyperparameter) == 1 && ~isempty(hyperparameter)
    imdl.hyperparameter.value = hyperparameter;
    % Remove function parameter because it will recurse
    try; imdl.hyperparameter = rmfield(imdl.hyperparameter,'func'); end
end


%% generate targets inside of rec_model
if isfield(imdl, 'rec_model')
    RecMdl = imdl.rec_model;   
else
    RecMdl = imdl.fwd_model;
end
MdlCtr = mean(RecMdl.nodes,1);

if isempty(xyzr_targets)
    % targets have not been defined, we'll generate them automatically 
    % by shrinking the model by the desired scaling factor and see the
    % descired number of targets inside the shrunk area
    %

    % first determine electrode level
    if isfield(RecMdl, 'mdl_slice_mapper') && isfield(RecMdl.mdl_slice_mapper, 'level')
        ElectrodeLevel = RecMdl.mdl_slice_mapper.level;    
    elseif isfield(RecMdl, 'mdl_slice_mapper') && isfield(RecMdl.mdl_slice_mapper, 'z_pts')
        ElectrodeLevel = RecMdl.mdl_slice_mapper.z_pts;    
        error('3D not supported yet, please seed targets manually via xyzr_targets');
    else
        if isfield(RecMdl.electrode(1), 'pos')
            elec_loc = cell2mat(cellfun(@(x) mean(x)', {RecMdl.electrode.pos}, 'uniformoutput', false))';
        elseif isfield(RecMdl.electrode(1), 'nodes')
            for i=1:length(RecMdl.electrode)
               enodesi = RecMdl.electrode(i).nodes;
               elec_loc(i,:) = mean( RecMdl.nodes( enodesi,:),1 );
            end  
        else
            error('not supported!');
        end
        ElectrodeLevel = mean(elec_loc,1);
    end
    
    % uniformly distribute N_DESIRED_TARGETS targets in the shrunken model 
    % at the level of the electrodes
    
    % first, shrink original model
    RecMdlShrunk = RecMdl;
    RecMdlShrunk.nodes = RecMdlShrunk.nodes - repmat(MdlCtr, size(RecMdlShrunk.nodes, 1), 1);
    RecMdlShrunk.nodes = RecMdlShrunk.nodes * ROI_SCALING_FACTOR;
    RecMdlShrunk.nodes = RecMdlShrunk.nodes + repmat(MdlCtr, size(RecMdlShrunk.nodes, 1), 1);

    % find boundary of shrunken model
    BoundaryNodes = find_boundary(RecMdlShrunk);
    % TODO: there's still a little problem here that we don't always
    % get a nicely oriented boundary, investigate this!
    Boundary = order_loop(RecMdlShrunk.nodes(BoundaryNodes, :));
    Boundary = [Boundary; Boundary(1,:)];
    % pp = fourier_fit(Boundary, 10); Boundary = fourier_fit(pp, linspace(0,1,20));

    % seed uniformly in rectangle and remove outliers
    AreaPoly = polyarea(Boundary(:,1), Boundary(:,2));
    Bounds = [max(Boundary); min(Boundary)];
    AreaRect = prod([Bounds(1,:) - Bounds(2,:)]);
    nUniformTgts = AreaRect * (N_DESIRED_TARGETS / AreaPoly);    

    % Size of the ROI in the two dimensions (x/y)...
    RoiSize = abs(diff(Bounds,[], 1));
    % ensure uniform spacing and scale according to differences in x/y size
    ScaleX = RoiSize(1) / RoiSize(2);
    ScaleY = 1./ScaleX;                 
    % create ROI centers
    [xx, yy] = meshgrid(linspace(Bounds(2,1), Bounds(1,1), ceil(ScaleX * sqrt(nUniformTgts))), ...
                        linspace(Bounds(2,2), Bounds(1,2), ceil(ScaleY * sqrt(nUniformTgts))));
    IsInside = inpolygon(xx(:), yy(:), Boundary(:,1), Boundary(:,2));
    nTgts = sum(IsInside);
    assert(abs((nTgts - N_DESIRED_TARGETS)/N_DESIRED_TARGETS) < 0.15, 'Cannot make desired number of targets');
    xx = xx(IsInside);
    yy = yy(IsInside);    
    %zz = ones(size(xx))*ElectrodeLevel(3);
else
    % targets specified from outside the function, assign them properly
    xx = xyzr_targets(1,:)';
    yy = xyzr_targets(2,:)';
    zz = xyzr_targets(3,:)';
end

% set target size relative to maximal model radius
BoundsFull = [max(RecMdl.nodes); min(RecMdl.nodes)];
Rmodel = (max(BoundsFull(1,:) - BoundsFull(2,:)))/2;    % maximal model radius
Rtarget = Rmodel * NORM_TGT_RADIUS;
rr = ones(size(xx))*Rtarget;

if ~isempty(xyzr_targets)
    try
        rr = xyzr_targets(4,:)';    % overwrite target radii if existing
    end
end
    
%% generate differential voltages for each conductivity target 

if elem_dim(imdl.fwd_model) == 3
    xyzr = [xx yy zz rr]';  
elseif elem_dim(imdl.fwd_model) == 2
    xyzr = [xx yy rr]';    
else
    error('unsupported dimensions');
end
try
    weighted_diff=imdl.inv_solve_complete_diff_GN_iter.weighted_diff;
catch
    disp('inv_solve_cole_diff_GN_iter: no weigthed difference is performed');
    weighted_diff=false;
end


scale = 0.1;
img = mk_image(imdl.fwd_model, scale);
img.c2f_scalar = 1.2*scale;
[vi1, xyzrOut1, c2f1] = my_simulate_movement(img, xyzr);

NotAssigned = ~ismember(xyzr', xyzrOut1', 'rows');
assert(sum(NotAssigned) == 0, 'Error: target(s) got missing...');


img = mk_image(imdl.fwd_model, 1.2*scale);
img.c2f_scalar = 2*scale;
[vi2, xyzrOut2, c2f2] = my_simulate_movement(img, xyzr);

NotAssigned = ~ismember(xyzr', xyzrOut2', 'rows');
assert(sum(NotAssigned) == 0, 'Error: target(s) got missing...');

if weighted_diff(1)
    alpha = (sum(vi1.*vi2,1)./sum(vi1.*vi1)); % factor calculation
else
    alpha = 1;
end
%% get reconstructed conductivity
try
    imdl.fwd_model = imdl.rec_model;
end

RM = get_fdRM(imdl,scale);
nb_elem = size(imdl.fwd_model.elems,1);
bkgnd = scale*ones(size(RM,1),1);
bkgnd(nb_elem+1:end,:) = 1.2*bkgnd(nb_elem+1:end,:);
img = mk_image(imdl.fwd_model, scale);
vh1 = fwd_solve(img);
vh1 = vh1.meas;
% vh1 = repmat(vh1,1,nTgts);

% img = mk_image(imdl.fwd_model, mean(alpha)*1.3*scale);
img.elem_data = alpha(1).*bkgnd(nb_elem+1:end,1);
% img.fwd_model.coarse2fine = alpha.*bkgnd(nb_elem+1:end,:);
% img.fwd_model = mdl_normalize(img.fwd_model,0);
% vh2 = calc_jacobian( img );
vh2 = fwd_solve(img);
vh2 = vh2.meas;
try
    with_v2=imdl.inv_solve_complete_diff_GN_iter.with_v2;
catch
    disp('inv_solve_complete_diff_GN_iter: no additional optimization regarding the measured voltages');
    with_v2=[false,0];
end
if with_v2(1)
    imgrs = -RM*[vh2-vh1-vi2+alpha(1).*vi1(:,1); with_v2(2)*(vh1-vi1(:,1)); with_v2(2)*(vh2-vi2(:,1))];
else
    imgrs = -RM*(vh2-vh1-vi2(:,1)+alpha(1).*vi1(:,1));
end
[step_all , step_half] = scale_to_fit_data(img, vh1, vh2,imgrs(:,1),bkgnd, with_v2,alpha(1),1);
imgrsNew = complex_truncation(bkgnd(:,1) + [repmat(step_all,nb_elem,1);repmat(step_all.*step_half,nb_elem,1)].*imgrs(:,1));

imgrs = calc_cond_change(imgrsNew,nb_elem);
%% generate individual target evaluation zones (TEZs)


imgrsNorm = imgrs;
imgrsNorm = imgrsNorm ./ repmat(max(imgrsNorm, [], 1), size(imgrsNorm,1), 1);
TEZs = double((imgrsNorm > TEZ_REGION_THRESH));

clear imgrsNorm;

% calculate volume/area of each element in RecMdl 
RecMdlVols = get_elem_volume(RecMdl);
ElemVols = spdiags(RecMdlVols, 0, length(RecMdlVols), length(RecMdlVols));
  
% TEZ matrices with proper normalization (as defined in publication)
%      z_i  =  a_i     *  c_i    (in publication)
Z = ElemVols * TEZs;
Vtez = sum(Z,1);                % volume/area of each TEZ
Vt = 4/3 * pi * rr.^3;          % volume of each target
% note: even for the 2D case we take the volume as a target double
% as thick than another but with same area has double the influence
K = Vtez ./ Vt';                % TEZ volume / target volume
% \tilde{z}_i = z_i * k          (in publication)
Z = Z./ repmat(sum(Z,1), size(TEZs,1), 1);
Ztilde = Z .* repmat(K, size(Z,1), 1);

% TEZ matrices kept for debugging / bw-compatility reasons
TEZsNonNorm = Z; 
% normalize ROIs 
TEZs = TEZs ./ repmat(sum(TEZs,1), size(TEZs,1), 1);  

%% get reconstructed conductivity ms in pure noise
nNoise = 20;
snr = linspace(1e3,5e3,nNoise);
Noise = zeros(1,nNoise);
for i = 1:nNoise
    vn1 = add_noise(snr(i),vi1(:,1),vh1);
    vn1 = vn1.meas;
    vn2 = add_noise(snr(i),vi2(:,1),vh2);
    vn2 = vn2.meas;
    if with_v2(1)
        imgrn = -RM*[vh2-vn2-vh1+alpha.*vn1;with_v2(2)*(vh1-vn1);with_v2(2)*(vh2-vn2)];
    else
        imgrn = -RM*(vh2-vn2-vh1+alpha.*vn1);
    end
    [step_all , step_half] = scale_to_fit_data(img, vh1, vh2,imgrn(:,1),bkgnd(:,1), with_v2,alpha(1),1);
    imgrnNew = complex_truncation(bkgnd + [repmat(step_all,nb_elem,1);repmat(step_all.*step_half,nb_elem,1)].*imgrn(:,1));
    imgrn = calc_cond_change(imgrnNew,nb_elem);
    
    VarPerPixel = diag(imgrn*imgrn');  
    Noise(:,i) = sqrt(Z' * VarPerPixel);
end
Noise = mean(Noise,2);

%% calculate target-wise distinguishability

Signal = diag(Ztilde(:,1)' * imgrs);    % numerator in equation (9)
%VarPerPixel = diag(RM*SigmaN*RM');  
%VarPerPixel_ = VarPerPixel(nb_elem+1:end)- VarPerPixel(1:nb_elem);
%Noise = sqrt(Z' * VarPerPixel_);      % denominator in equation (9)  

% TODO: make sure how to get difference conductivities of noise
%Noise = sqrt(Z' * VarPerPixel(nb_elem+1:end))-sqrt(Z' * VarPerPixel(1:nb_elem));

% get the average of all SNRs
SNRs = Signal ./ Noise; 
if find_min
    SNRmean = -mean(SNRs); 
    SNR_ = -SNRmean;
else
    SNRmean = mean(SNRs); 
    SNR_ = SNRmean;
end
SE = std(SNRs);

if isnan(SNRmean)
    SNRmean = -inf;
elseif isinf(SNRmean)
    keyboard;
end

try
  eidors_msg('SNR = %e (hp=%e)', SNR_, imdl.hyperparameter.value, 1);
end


%% debug info: fraction of amplitude response (AR) inside each TEZ
if nargout >= 3
    % target-wise amplitude response inside TEZ: in terms of energy
    ArInFrac = diag((TEZsNonNorm'*abs(RM*vd).^2) ./ (repmat(RecMdlVols, 1, size(TEZs,2))'*abs(RM*vd).^2));    
    debug.ArInFrac = ArInFrac;
    
    debug.SNRs = SNRs;
    debug.Signal = Signal;
    debug.Noise = Noise;
    debug.TEZs = TEZs;
    debug.TEZsNoNorm = TEZsNonNorm;
    
    debug.RoiBounds = Bounds;
    debug.RoiBoundary = Boundary;
    debug.MdlCtr = MdlCtr;
    debug.BoundsFull = BoundsFull;
    
    imgrsNorm = imgrs;
    imgrsNorm = imgrsNorm ./ repmat(max(imgrsNorm, [], 1), size(imgrsNorm,1), 1);
    debug.meanInNormImg = diag(TEZs'*imgrsNorm);
end


%% visualize for debug purposes
if doPlot
    fig = figure;
    set(fig, 'position', [  182         700        1531         313]);
    
    if isfield(imdl.fwd_model, 'coarse2fine')
        MapTgts2Img = imdl.fwd_model.coarse2fine'*c2f;
    else
        try
            imdl.fwd_model.coarse2fine = mk_coarse_fine_mapping(imdl.fwd_model, imdl.rec_model);
            MapTgts2Img = imdl.fwd_model.coarse2fine'*c2f;
        catch
            warning('Unable to make c2f mapping');
            MapTgts2Img = c2f;
        end
    end
    MapTgts2Img = MapTgts2Img ./ repmat(sum(MapTgts2Img,2), 1, size(MapTgts2Img, 2));
    MapTgts2Img(isnan(MapTgts2Img(:))) = 0;
    
    img = mk_image(RecMdl, nan);
    img.elem_data = MapTgts2Img * SNRs;
    SnrImg = calc_slices(img);
    
    img.elem_data = nan(size(img.elem_data));
    img.elem_data = MapTgts2Img * Signal;
    AmpSens = img.elem_data;
    AmpSensImg = calc_slices(img);
    
    img.elem_data = nan(size(img.elem_data));
    img.elem_data = MapTgts2Img * Noise;
    NoiseSens = img.elem_data;
    NoiseSensImg = calc_slices(img);
    
    sp1 = subplot(131);
    imagescnan(SnrImg);
    title(['SNR image: ', num2str(SNRmean, '%0.2d')]);
    colorbar; colormap jet;

    % plot signal sensitivity image
    sp2 = subplot(132);
    imagescnan(AmpSensImg);
    title(['Amplitude response: ', num2str(nanmean(AmpSens(:)), '%0.2d')]);
    colorbar;

    % plot noise sensitivity image
    sp3 = subplot(133);
    imagescnan(NoiseSensImg);
    title(['Noise sensitivity: ', num2str(nanmean(NoiseSens(:)), '%0.2d')]);
    colorbar;

    if doPlot < 2
        linkaxes([sp1 sp2 sp3], 'xy');
    end
   
    if doPlot == 2
       % some extra debugging 
       sp1 = subplot(131);
       imgTmp = img;
       iTgt = round(size(vd,2)/2);
       imgTmp.elem_data = RM*vd(:,iTgt);
       hh = show_fem(imgTmp);
       set(hh, 'edgecolor', 'none');
       hold on; 
       circle(xyzr(1:2, iTgt), xyzr(3, iTgt), 100, 'k');
       circle(xyzr(1:2, iTgt), Rtarget, 100, 'k');
       axis equal;
    end
    
    
    if doPlot == 3
       sp3 = subplot(133); 
       cla;
       hh = show_fem(RecMdl);
       set(hh, 'edgecolor', 'none')
       hold on;
       plot(Boundary(:,1), Boundary(:,2), '.-k');
       plot(xyzr(1,:), xyzr(2,:), '.r');
       plot(xyzr(1,:), xyzr(2,:), 'or');
    end
    
    if doPlot == 4
        % show pixel-wise noise sensitivity
       sp2 = subplot(132);
       cla;
       imgN = mk_image(imdl,0);
       imgN.elem_data = VarPerPixel;
       imgN.elem_data(sum(MapTgts2Img,2) == 0) = 0;
       imagescnan(calc_slices(imgN));
       title(['Pixel-wise NoiseSens: ', num2str(nanmean(imgN.elem_data(:)), '%0.2d')]);
       colorbar;        
    end
    
    if doPlot == 5
        sp1 = subplot(131);
        hist(SNRs); 
        xlim([0 max(xlim())]);
    end
    
end

end


% function [step_all,step_half] = scale_to_fit_data(img, inv_model,data1, data2,delta_sol,sol, with_v2, alpha)
%    
%       eidors_msg('inv_solve_complete_diff_GN_iter: Calculating optimal step size to fit data',2);
%       % options for fmincon
%       opts = optimoptions(@fmincon,'Display','iter');
% %         opts = optimoptions('ga','Display','iter');
%     
%       % range for fminbnd
%         try
%         range = inv_model.inv_solve_complete_diff_GN_iter.bounds;
%         catch
%         range = [1e-6 1e3];
%         end
% 
%         aux=-real(sol)./real(delta_sol);
%         x_range2=min(aux(aux>range(1)));
% 
%         if isempty(x_range2) 
%         x_range2=range(2);
%         elseif x_range2>range(2)
%         x_range2=range(2);
%         end
% 
%         clearvars x;
%         noncon = @(x)calc_cond_bound(sol,delta_sol,x);
%         step = fmincon(@(x) to_optimize_all1(img,data1,data2,delta_sol,sol, x,with_v2, alpha), ...
%                    [1,1],[],[],[],[], range(1).*[1,1],x_range2.*[1,1],noncon, opts);
%                         
%         step_all = step(1);
%         step_half = step(2);
% 
% end   
function [step_all,step_half] = scale_to_fit_data(img,data1, data2,delta_sol,sol, with_v2, alpha,nTgts)
    range = [1e-8 4];
    step_all  = zeros(1,nTgts);
    step_half = zeros(1,nTgts);
    opt = optimset('Display','iter');
    for i = 1:size(delta_sol,2)
        clearvars x;
        [step_all(i),fval1] = fminbnd(@(x) to_optimize_all1(img,data1(:,i),data2(:,i),delta_sol(:,i),sol(:,i), x,with_v2, alpha(i)), ...
        range(1), -0.1/min(delta_sol(:,i)),opt);

        clearvars x;
        [step_half(i),fval2] = fminbnd(@(x) to_optimize_half(img,data1(:,i),data2(:,i),step_all(i)*delta_sol(:,i),sol(:,i), x,with_v2, alpha(i)), ...
        range(1), range(2),opt);
        % If the second step-size line-search didn't improve anything, then
        % ignore it
        if fval2>fval1
            step_half(i)=1;
        end
        disp(i);
    end   
end

function out = to_optimize_all(img, data1, data2,delta_sol,sol, x,with_v2, alpha)
    %Function f_global in the Masterarbeit. The subtleties (mixed 
    %minimization and normalization are considered)
    nb_elem=size(sol,1)/2;
    img.elem_data= - complex_truncation(sol(1:nb_elem)+x*delta_sol(1:nb_elem)); %rajout de complex truncation pour eviter les conductivités negatives
    Psi1=fwd_solve(img);
    Psi1=Psi1.meas;

%     img.elem_data=complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x*delta_sol((nb_elem+1):2*nb_elem,:));
%     Psi2=fwd_solve(img);
%     Psi2=Psi2.meas;
%     img.elem_data=alpha.*complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x*delta_sol((nb_elem+1):2*nb_elem,:));
%     Psi2w=fwd_solve(img);
%     Psi2w=Psi2w.meas;
    
    % Second term due to the mixed minimization
%     if with_v2(1)
%             corr=with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
%     else
%         corr=0;
%     end
%     out=norm(Psi2w-(Psi1)-data2+(alpha.*data1))+corr;

    out = norm((Psi1-data1));
    if isnan(out)
        out=x*10e6;
    end
end
function out = to_optimize_all1(img, data1, data2,delta_sol,sol, x,with_v2, alpha)
    %Function f_global in the Masterarbeit. The subtleties (mixed 
    %minimization and normalization are considered)
    nb_elem=size(sol,1)/2;
    img.elem_data=complex_truncation(sol(1:nb_elem)+x*delta_sol(1:nb_elem)); %rajout de complex truncation pour eviter les conductivités negatives
    Psi1=fwd_solve(img);
    Psi1=Psi1.meas;

    img.elem_data=complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x*delta_sol((nb_elem+1):2*nb_elem,:));
    Psi2=fwd_solve(img);
    Psi2=Psi2.meas;
    img.elem_data=alpha.*complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x*delta_sol((nb_elem+1):2*nb_elem,:));
    Psi2w=fwd_solve(img);
    Psi2w=Psi2w.meas;
%     cond1 = complex_truncation(sol(1:nb_elem)+x(1)*delta_sol(1:nb_elem));
%     cond2 = alpha.*complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x(2)*delta_sol((nb_elem+1):2*nb_elem,:));
%     img.elem_data = cond2-cond1;
%     dPsiw = fwd_solve(img);
%     dPsiw = dPsiw.meas;
%     Second term due to the mixed minimization
    if with_v2(1)
            corr=with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
    else
        corr=0;
    end
    out=norm(Psi2w-Psi1-data2+(alpha.*data1))+corr;

%     out = norm((Psi1-data1));
    if isnan(out)
        out=x*10e6;
    end
end    
function out = to_optimize_half(img, data1, data2,delta_sol,sol,x,with_v2,alpha)
%Function f_half in the Masterarbeit. The subtleties (mixed 
%minimization and normalization are considered)
    nb_elem=size(sol,1)/2;
    img.elem_data=complex_truncation(sol(1:nb_elem,:)+delta_sol(1:nb_elem,:));
    Psi1=fwd_solve(img);
    Psi1=Psi1.meas;

    img.elem_data=complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x*delta_sol((nb_elem+1):2*nb_elem,:));
    Psi2=fwd_solve(img);
    Psi2=Psi2.meas;
    img.elem_data=alpha.*complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x*delta_sol((nb_elem+1):2*nb_elem,:));
    Psi2w=fwd_solve(img);
    Psi2w=Psi2w.meas;

    if with_v2(1)
        corr=with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
    else
        corr=0;
    end

    out=norm(Psi2w-Psi1-data2+(alpha.*data1))+corr;        

    if isnan(out)
        out=x*10e6;
    end
end
   
function cond_change=calc_cond_change(sol,nb_elem)
% Function giving the conductivity change between gamma^(high) and gamma^(low)
cond_change=sol(nb_elem+1:2*nb_elem,:)-sol(1:nb_elem,:);
end

function [c,ceq] =  calc_cond_bound(sol,delta_sol,x)

nb_elem = size(sol,1)/2;
c = -real(sol+[x(1)*delta_sol(1:nb_elem);x(1)*x(2)*delta_sol(nb_elem+1:end)]);
ceq = [];
end
