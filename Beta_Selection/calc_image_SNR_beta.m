function [SNRmean,imgrs, SE, debug] = calc_image_SNR_beta(imdl, hp1,hp2, doPlot)
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
% get conductivities of tissues 
freq_1 = 100e3; freq_2 = 700e3;
tissues_1 = getTissues_mod(freq_1);
tissues_2 = getTissues_mod(freq_2); 
s_target_1 = round(tissues_1.Muscle_combined,4,'decimal');
s_target_2 = round(tissues_2.Muscle_combined,4,'decimal'); 
s_bkgnd_1 = round(tissues_1.LungDeflated_combined,4,'decimal'); 
s_bkgnd_2 = round(tissues_2.LungDeflated_combined,4,'decimal');
% s_bkgnd_1 = round(tissues_1.Muscle_combined,4,'decimal');
% s_bkgnd_2 = round(tissues_2.Muscle_combined,4,'decimal'); 
% s_target_1 = round(tissues_1.LungDeflated_combined,4,'decimal'); 
% s_target_2 = round(tissues_2.LungDeflated_combined,4,'decimal');
% Specify reference homogenous conductivities and conductivities for 
% targets 
% Remark: The conductivity for low and high frequency should be complex,
% otherwise the background difference is just a real value, this does not
% fullfill the requierements for fdEIT.
img = mk_image(imdl.fwd_model, s_bkgnd_1);
img.c2f_scalar = s_target_1;
[vi1, xyzrOut1, c2f1] = my_simulate_movement(img, xyzr);

NotAssigned = ~ismember(xyzr', xyzrOut1', 'rows');
assert(sum(NotAssigned) == 0, 'Error: target(s) got missing...');


img = mk_image(imdl.fwd_model, s_bkgnd_2);
img.c2f_scalar = s_target_2;
[vi2, xyzrOut2, c2f2] = my_simulate_movement(img, xyzr);

NotAssigned = ~ismember(xyzr', xyzrOut2', 'rows');
assert(sum(NotAssigned) == 0, 'Error: target(s) got missing...');

if weighted_diff(1)
    alpha = (sum(vi1.*vi2,1)./sum(vi1.*vi1)); % factor calculation
else
    alpha = 1;
end
try 
    mode = imdl.inv_solve_complete_diff_GN_iter.mode;
catch
    mode=1;
end
try
    with_v2=imdl.inv_solve_complete_diff_GN_iter.with_v2;
    if size(with_v2,2) == 1&& with_v2==0
         with_v2 = [false,0];
    elseif size(with_v2,2) == 1&& with_v2~=0
         with_v2 = [true,with_v2];
    end
catch
    disp('inv_solve_complete_diff_GN_iter: no additional optimization regarding the measured voltages');
    with_v2=[false,1];
end
%% get reconstructed conductivity
try
    imdl.fwd_model = imdl.rec_model;
end

scale=1;
imdl.hyperparameter.value = hp1;
imdl.inv_solve_complete_diff_GN_iter.mode=1;
RM1 = get_fdRM(imdl,scale*s_bkgnd_1,scale*s_bkgnd_2,mean(alpha));
imdl.hyperparameter.value = hp2;
imdl.inv_solve_complete_diff_GN_iter.mode=2;
RM2 = get_fdRM(imdl,scale*s_bkgnd_1,scale*s_bkgnd_2,mean(alpha));
nb_elem = size(imdl.fwd_model.elems,1);
bkgnd = scale*s_bkgnd_1*ones(size(RM1,1),nTgts);
bkgnd(nb_elem+1:end,:) = scale*s_bkgnd_2;
% Remark: the "homogeneous" voltages are not the same for all targets,
% we perform fdEIT and also the low conductivity measurement changes for
% each target.

% voltage measurement corresponding to Psi1 and Psi2w in algorithm
img = mk_image(imdl.fwd_model, scale*s_bkgnd_1);
vh1 = fwd_solve(img);
vh1 = vh1.meas;
vh1 = repmat(vh1,1,nTgts);

img = mk_image(imdl.fwd_model,scale*mean(alpha).*s_bkgnd_2);
vh2 = fwd_solve(img);
vh2 = vh2.meas;
vh2 = repmat(vh2,1,nTgts);

img = mk_image(imdl.fwd_model,scale*s_bkgnd_2);
vh2_ = fwd_solve(img);
vh2_ = vh2_.meas;
vh2_ = repmat(vh2_,1,nTgts);
% clear vh2;
% for xx = 1:length(vh)
%     vh2(:,xx) = vh(xx).meas;
% end
% weighted difference
imgrs1 = -RM1*(vh2-vh1-vi2+alpha.*vi1);
% absolute
imgrs2 = -RM2*[(vh1-vi1); (vh2_-vi2)];


% Remark: The step size calculation takes a long time, but we assume the 
% step size does not differ much in between targets (targets are relatively
% small).
% Therefore, we can use the background conductivity of f1 and f2 to
% calculate a single step size for all targets. -> great performance
% increase
% step size
[step_all1 , step_half1] = scale_to_fit_data(img, vh1, vh2,imgrs1,bkgnd, with_v2,alpha,nTgts,10,1);
[step_all2 , step_half2] = scale_to_fit_data(img, vh1, vh2,imgrs2,bkgnd, with_v2,alpha,nTgts,10,2);

%  both weighted different and absolute rconstruction
imgrsNew1 = [repmat(step_all1,nb_elem,1);repmat(step_all1.*step_half1,nb_elem,1)].*imgrs1;
imgrs1 = calc_cond_change(imgrsNew1,nb_elem);
imgrsNew2 = [repmat(step_all2,nb_elem,1);repmat(step_all2.*step_half2,nb_elem,1)].*imgrs2;
imgrs2 = calc_cond_change(imgrsNew2,nb_elem);
% imgrs = calc_cond_change(imgrs,nb_elem);
%% get reconstructed conductivity ms in pure noise
if isfield(imdl,'beta')
    beta = imdl.beta;
else
    beta = 0:0.1:1;
end
nNoise = 20;
snr = 50;
Noise = zeros(nTgts,nNoise);
step_alln1 = step_all1;
step_halfn1 = step_half1;
step_alln2 = step_all2;
step_halfn2 = step_half2;
for b = 1:length(beta)
    imgrs = imgrs1+beta(b)*imgrs2;
%% generate individual target evaluation zones (TEZs)    
imgrsNorm = imgrs;
imgrsNorm = imgrsNorm ./ repmat(max(real(imgrsNorm), [], 1), size(imgrsNorm,1), 1);
% imgrsNorm = imgrsNorm ./ repmat(max(real(imgrsNorm), [], 1), size(imgrsNorm,1), 1);
TEZs = double((real(imgrsNorm) > TEZ_REGION_THRESH));

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

for i = 1:nNoise
    noise = randn(2*size(vi1,1),2*nTgts);
    SNR = norm(vi1-vh1)/norm(noise)/snr;
    vn1 = vi1+SNR*noise(1:size(vi1,1),1:nTgts);
    vn2 = vi2+SNR*noise(size(vi1,1)+1:end,nTgts+1:end);
%     if weighted_diff(1)
%         alpha = (sum(vn1.*vn2,1)./sum(vn1.*vn1)); % factor calculation
%     else
%         alpha = 1;
%     end
  
    imgrn1 = -RM1*(vh2-vh1-vn2+alpha.*vn1);
    imgrn2 = -RM2*[ (vh1-vn1); (vh2_-vn2)];
%     if mod(i,5)== 0
%     [step_alln , step_halfn] = scale_to_fit_data(img, vh1, vh2,imgrn,bkgnd, with_v2,alpha,nTgts,1,mode);
%     end
    imgrnNew1 =  [repmat(step_alln1,nb_elem,1);repmat(step_alln1.*step_halfn1,nb_elem,1)].*imgrn1;
    imgrn1 = calc_cond_change(imgrnNew1,nb_elem);
    imgrnNew2 =  [repmat(step_alln2,nb_elem,1);repmat(step_alln2.*step_halfn2,nb_elem,1)].*imgrn2;
    imgrn2 = calc_cond_change(imgrnNew2,nb_elem);
    imgrn = imgrn1+beta(b)*imgrn2-imgrs;
    VarPerPixel = diag(real(imgrn)*real(imgrn)');  
    Noise(:,i) = Z' * VarPerPixel;
    clear VarPerPixel imgrn imgrnNew vn1 vn2
end
%% calculate target-wise distinguishability

Signal = diag(Ztilde' * real(imgrs));    % numerator in equation (9)
Noise = sqrt(sum(Noise,2)/(nNoise-1));

% get the average of all SNRs
SNRmean(b) = mean(Signal ./ Noise);
end
 

if isnan(SNRmean)
    SNRmean = -inf;
elseif isinf(SNRmean)
    keyboard;
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
function [step_all,step_half] = scale_to_fit_data(img,data1, data2,delta_sol,sol, with_v2, alpha,nTgts,nSteps,mode)
    range = [1e-8 4];
    opt = optimset('Display','iter','MaxIter',30);
    nb_elem = size(sol,1)/2;
    blank = idivide(nTgts,int32(nSteps));
    a = 1;
    for i =1:blank:nTgts-blank
    bound1 =-min(real(sol(:,i)))/min(real(delta_sol(:,i)));
    clearvars x;
    
    [step_all(a),fval1] = fminbnd(@(x) to_optimize_all(img,data1(:,i),data2(:,i),delta_sol(:,i),sol(:,i), x,0,with_v2, alpha(i),mode), ...
    range(1),bound1);

    clearvars x;
    bound2 =-min(real(sol(nb_elem+1:end,i)))/min(real(step_all(a)*delta_sol(nb_elem+1:end,i)));
    [step_half(a),fval2] = fminbnd(@(x) to_optimize_half(img,data1(:,i),data2(:,i),step_all(a)*delta_sol(:,i),sol(:,i), x,0,with_v2, alpha(i),mode), ...
    range(1), bound2);
    % If the second step-size line-search didn't improve anything, then
    % ignore it
%     if fval2>fval1
%         step_half(a)=1;
%     end
%     step_half(a) = 1;
    a = a+1;
    end
    step_all = kron(step_all,ones(1,blank));
    step_half = kron(step_half,ones(1,blank));
    if size(step_all,2)~=nTgts
        step_all(size(step_all,2):nTgts)=step_all(size(step_all,2));
        step_half(size(step_half,2):nTgts)=step_half(size(step_half,2));
    end
%     step_half = repmat(step_half,1,nTgts);
    
end

function out = to_optimize_all(img, data1, data2,delta_sol,sol, x,normalize,with_v2, alpha,mode)
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
    
    % Second term due to the mixed minimization

    % Normalized minimization
    
    if normalize
        if mode == 1
            out = norm(Psi2w./(Psi1)-data2./(alpha.*data1))+with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
        elseif mode == 2
%             out = (norm((Psi1-data1))+norm((Psi2-data2)));
            out = (norm((Psi1-data1)));
%             out = norm(Psi2w./(Psi1)-data2./(alpha.*data1))+0.2*(norm((Psi1-data1))+norm((Psi2-data2)));
        elseif mode == 3
            out = norm(Psi2w./(Psi1)-data2./(alpha.*data1));
        end
    else
        if mode == 1
            out = norm(Psi2w-Psi1-data2+alpha.*data1)+with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
        elseif mode == 2
            out = norm((Psi1-data1));
%             out = norm(Psi2w-Psi1-data2+alpha.*data1);
        elseif mode == 3
            out = norm(Psi2w+Psi1-data2-alpha.*data1);
        end
    end
  
    if isnan(out)
        out=x*10e6;
    end
end
    
function out = to_optimize_half(img, data1, data2,delta_sol,sol,x,normalize,with_v2,alpha,mode)
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
    
    if normalize
        if mode == 1
            out = norm(Psi2w./(Psi1)-data2./(alpha.*data1))+with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
        elseif mode == 2
            out = (norm((Psi1-data1))+norm((Psi2-data2)));
            %out = norm(Psi2w./(Psi1)-data2./(alpha.*data1))+(norm((Psi1-data1))+norm((Psi2-data2)));
        elseif mode == 3
            out = norm(Psi2w./(Psi1)-data2./(alpha.*data1));
        end
    else
        if mode == 1
            out = norm(Psi2w-Psi1-data2+alpha.*data1)+with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
        elseif mode == 2
            out = (norm((Psi1-data1))+norm((Psi2-data2)));
%             out = norm(Psi2w-Psi1-data2+alpha.*data1);
        elseif mode == 3
            out = norm(Psi2w+Psi1-data2-alpha.*data1);
        end
    end
  
    
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
function [step_all,step_half] = scale_to_fit_data1(img,data1, data2,delta_sol,sol, with_v2, alpha,nTgts,nSteps)
    range = [1e-8 4];
    opt = optimset('Display','iter','MaxIter',30);
    nb_elem = size(sol,1)/2;
    blank = idivide(nTgts,int32(nSteps));
    a = 1;
    for i =1:blank:nTgts-1
    bound1 = max(real(sol(:,i)))/max(real(delta_sol(:,i)));
    clearvars x;
    
    [step_all(a),fval1] = fminbnd(@(x) to_optimize_all11(img,data1(:,i),data2(:,i),delta_sol(:,i),sol(:,i), x,with_v2, alpha(i)), ...
    range(1),bound1);

    clearvars x;
    bound2 =max(real(sol(nb_elem+1:end,i)))/max(real(delta_sol(nb_elem+1:end,i)));
    [step_half(a),fval2] = fminbnd(@(x) to_optimize_half11(img,data1(:,i),data2(:,i),delta_sol(:,i),sol(:,i), x,with_v2, alpha(i)), ...
    range(1), bound2);
    % If the second step-size line-search didn't improve anything, then
    % ignore it
%     if fval2>fval1
%         step_half(a)=1;
%     end
%     step_half(a) = 1;
    a = a+1;
    end
    step_all = kron(step_all,ones(1,blank));
    step_half = kron(step_half,ones(1,blank));
    if size(step_all,2)~=nTgts
        step_all(size(step_all,2):nTgts)=step_all(size(step_all,2));
        step_half(size(step_half,2):nTgts)=step_half(size(step_half,2));
    end
%     step_half = repmat(step_half,1,nTgts);
    
end

function out = to_optimize_all11(img, data1, data2,delta_sol,sol, x,with_v2, alpha)
    %Function f_global in the Masterarbeit. The subtleties (mixed 
    %minimization and normalization are considered)
    nb_elem=size(sol,1)/2;
    img.elem_data= complex_truncation(sol(1:nb_elem)+x*delta_sol(1:nb_elem)); %rajout de complex truncation pour eviter les conductivités negatives
    Psi1=fwd_solve(img);
    Psi1=Psi1.meas;

    out = norm((Psi1-data1));
    if isnan(out)
        out=x*10e6;
    end
end
function out = to_optimize_half11(img, data1, data2,delta_sol,sol,x,with_v2,alpha)
%Function f_half in the Masterarbeit. The subtleties (mixed 
%minimization and normalization are considered)
    nb_elem=size(sol,1)/2;
   
    img.elem_data=complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x*delta_sol((nb_elem+1):2*nb_elem,:));
    Psi2=fwd_solve(img);
    Psi2=Psi2.meas;


    out=norm(Psi2-data2);        

    if isnan(out)
        out=x*10e6;
    end
end