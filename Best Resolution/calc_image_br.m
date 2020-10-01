function BR = calc_image_br(imdl,hyperparameter)

if nargin>=2 && numel(hyperparameter) == 1 && ~isempty(hyperparameter)
    imdl.hyperparameter.value = hyperparameter;
end
try
    nb_noise = imdl.nb_noise;
    SNR = imdl.NoiseSNR;
catch
    nb_noise = 50;
    SNR = linspace(50,60,nb_noise);
end
try
    with_v2=imdl.inv_solve_complete_diff_GN_iter.with_v2;
catch
    disp('inv_solve_complete_diff_GN_iter: no additional optimization regarding the measured voltages');
    with_v2=[false,0];
end
try
    weight = imdl.inv_solve_complete_diff_GN_iter.weighted_diff;
catch
    weight = false;
end
try
    mode = imdl.inv_solve_complete_diff_GN_iter.mode;
catch
    mode = 1;
end
FwdMdlVols = get_elem_volume(imdl.fwd_model);
RecMdlVols = get_elem_volume(imdl.rec_model);
radius = sqrt(max(FwdMdlVols)/pi);
% select_fcn = inline(strcat('(x+0.5).^2+(y-0).^2<',num2str(radius),'^2'),'x','y','z');
select_fcn = inline(strcat('(x-0.5).^2+(y).^2<',num2str(radius),'^2'),'x','y','z');
% select_fcn = inline(strcat('(x-0.5).^2+(y).^2<0.1^2'),'x','y','z');
% get conductivities of tissues 
freq_1 = 100e3; freq_2 = 700e3;
tissues_1 = getTissues_mod(freq_1);
tissues_2 = getTissues_mod(freq_2); 
s_target_1 = round(tissues_1.Muscle_combined,4,'decimal');
s_target_2 = round(tissues_2.Muscle_combined,4,'decimal'); 
s_bkgnd_1 = round(tissues_1.LungDeflated_combined,4,'decimal'); 
s_bkgnd_2 = round(tissues_2.LungDeflated_combined,4,'decimal');

img = mk_image(imdl.fwd_model,s_bkgnd_1);
vr1 = fwd_solve(img);
vr1 = vr1.meas;
Elem = elem_select(img.fwd_model, select_fcn);
% img.elem_data =  s_bkgnd_1 ;
% img.elem_data(Elem~= 0) =  s_target_1;
img.elem_data(Elem==max(Elem)) =  s_target_1;
% img.elem_data(img.elem_data~=max(img.elem_data)) = 1;
v1 = fwd_solve(img);
v1 = v1.meas;

img = mk_image(imdl.fwd_model,s_bkgnd_2);
vr2 = fwd_solve(img);
vr2 = vr2.meas;
% img.elem_data =  s_bkgnd_2;
% img.elem_data(Elem ~= 0) =  s_target_2;
img.elem_data(Elem==max(Elem)) =  s_target_2;
% img.elem_data(img.elem_data~=max(img.elem_data)) = 1.2;
v2 = fwd_solve(img);
v2 = v2.meas;

A_0 = sum(FwdMdlVols);


clear vn1 vn2;
for i = 1:nb_noise
    n1 = add_noise(SNR(i),v1,vr1);
    n2 = add_noise(SNR(i),v2,vr2);
    vn1(:,i) = n1.meas;
    vn2(:,i) = n2.meas;
end

if weight
    alpha = diag(v1'*v2)./diag(v1'*v1);
else
    alpha = ones(nb_noise,1);
end
scale = 1;
img = mk_image(imdl.fwd_model,scale*s_bkgnd_1);
vr1 = fwd_solve(img);
vr1 = vr1.meas;

img = mk_image(imdl.fwd_model,scale*s_bkgnd_2);
vr2_ = fwd_solve(img);
vr2_ = vr2_.meas;

img = mk_image(imdl.fwd_model,scale*s_bkgnd_2);
img.elem_data = img.elem_data.*mean(alpha);
vr = fwd_solve(img);
vr2 = zeros(size(vr1,1),nb_noise);
for i = 1:nb_noise
    vr2(:,i) = vr.meas;
end

rmdl = imdl;
rmdl.fwd_model = imdl.rec_model;
nb_elem = size(rmdl.fwd_model.elems,1);
% init = 0.1+0.01i;
% step_optim=fminbnd(@(x) to_optimize_half(img, vn1, vn2,(0.1+0.01i)*ones(2*nb_elem,1),zeros(nb_elem,1), x,0,with_v2,1), ...
% 1e-5, 10);

img =  mk_image(rmdl.fwd_model,s_bkgnd_1);
RM = get_fdRM(rmdl,scale*s_bkgnd_1,scale*s_bkgnd_2,mean(alpha));
% RM = get_fdRM(rmdl,init,step_optim*init);
bkgnd = scale*s_bkgnd_1*ones(size(RM,1),nb_noise);
bkgnd(nb_elem+1:end,:) = scale*s_bkgnd_2*ones(size(RM,1)/2,nb_noise);


if mode==1
    imgrs = -RM*[vr2-vr1-vn2+mean(alpha).*vn1; with_v2(2)*(vr1-vn1); with_v2(2)*(vr2_-vn2)];
elseif mode==2
    imgrs = -RM*[(vr1-vn1); (vr2_-vn2)];
end
% Step = scale_to_fit_data(img, imdl,vn1, vn2,imgrs,bkgnd, with_v2, alpha);
% imgrsNew = stepsize_optimization(imdl,v1,v2,sol,delta_sol,alpha,'original');
% imgrsNew = calc_cond_change([repmat(Step(1,:),nb_elem,1);repmat(Step(2,:),nb_elem,1)].*imgrs,nb_elem); 
% imgrsNew = calc_cond_change(Step.*imgrs,nb_elem); 
imgrsNorm = calc_cond_change(imgrs,nb_elem);
% img_diff = imgrsNew - (sigma_freq2-sigma_freq1);
% imgrsNorm = imgrsNew - min(imgrsNew,[],1);
% imgrsNorm = imgrs - repmat(mean(imgrs), size(imgrs, 1), 1);
imgrsNorm = imgrsNorm -min(real(imgrsNorm));
imgrsNorm = imgrsNorm ./ repmat(max(real(imgrsNorm), [], 1), size(imgrsNorm,1), 1);
% imgrsNorm = normalize(real(imgrs),'range');
SelectVols = repmat(RecMdlVols,1,nb_noise).*(imgrsNorm>0.25);
clear imgrsNorm;
A_z = mean(sum(SelectVols));
BR = sqrt(A_z/A_0);

try
  eidors_msg('Blur Radius = %e (hp=%e)', BR, imdl.hyperparameter.value, 1);
end
end



function Step = scale_to_fit_data(img, inv_model,data1, data2,delta_sol,sol, with_v2, alpha)
    % Function to find the step size to multiply sol by to best fit data
    Step = zeros(2,size(delta_sol,2));
    % By default, the step size is computed.
    do_step   = true;
    try 
        do_step = inv_model.inv_solve_complete_diff_GN_iter.calc_step_size; 
    end

    if do_step
    %if the step-size must be calculated
        eidors_msg('inv_solve_complete_diff_GN_iter: Calculating optimal step size to fit data',2);
        % options for fminbnd
        try 
            opt = inv_model.inv_solve_complete_diff_GN_iter.fminbnd;
        catch
            opt.Display = 'iter';
            opt.ConstraintTolerance = 1e-10;
        end
        % range for fminbnd
        try
            range = inv_model.inv_solve_complete_diff_GN_iter.bounds;
        catch
            range = [1e-6 1e3];
        end

        try
            normalize=inv_model.inv_solve_complete_diff_GN_iter.normalize;
        catch
            disp('normalize=0');
            normalize=0;
        end
         hp = calc_hyperparameter( inv_model );
        for i = 1:size(delta_sol,2)
            clearvars x;
%             noncon = @(x)calc_cond_bound(sol(:,i),delta_sol(:,i),x);
%             step1 = -min(sol(:,i))/min(delta_sol(:,i));
%             step2 = 0.2*min(sol(:,i))/max(delta_sol(1:size(delta_sol,1)/2,i));
%             Step(1,i) = min([step1,step2]);
             bound1 =-min(real(sol(:,i)))/min(real(delta_sol(:,i)));
            Step(1,i) = fminbnd(@(x) to_optimize_all(img,data1(:,i),data2(:,i),delta_sol(:,i),sol(:,i), x,normalize,with_v2, alpha(i)), ...
                    range(1), bound1, opt);
            bound2 =-min(real(sol(nb_elem+1:end,i)))/min(real(Step(1,i)*delta_sol(nb_elem+1:end,i)));
            Step(2,i) = fminbnd(@(x) to_optimize_half(img,data1(:,i),data2(:,i),Step(1,i)*delta_sol(:,i),sol(:,i), x,normalize,with_v2, alpha(i)), ...
                range(1), bound2, opt);
            disp(strcat('Current hp:',num2str(hp)));
            disp(strcat('Process:',num2str(i),'/',num2str(size(delta_sol,2))));
        end
    end           
end
function out = to_optimize_step(img, data1, data2,delta_sol,sol, x,normalize,with_v2, alpha)
    %Function f_global in the Masterarbeit. The subtleties (mixed 
    %minimization and normalization are considered)
    nb_elem=size(sol,1)/2;
    img.elem_data=complex_truncation(sol(1:nb_elem)+x(1)*delta_sol(1:nb_elem)); %rajout de complex truncation pour eviter les conductivités negatives
    Psi1=fwd_solve(img);
    Psi1=Psi1.meas;

    img.elem_data=complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x(1)*x(2)*delta_sol((nb_elem+1):2*nb_elem,:));
    Psi2=fwd_solve(img);
    Psi2=Psi2.meas;
    img.elem_data=alpha.*complex_truncation(sol((nb_elem+1):2*nb_elem,:)+x(1)*x(2)*delta_sol((nb_elem+1):2*nb_elem,:));
    Psi2w=fwd_solve(img);
    Psi2w=Psi2w.meas;
    
    % Second term due to the mixed minimization
    if with_v2(1)
            corr=with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
    else
        corr=0;
    end

    % Normalized minimization
    if normalize
        out=norm(Psi2w./(Psi1)-data2./(alpha.*data1))+corr;
    else
        out=norm(Psi2w-(Psi1)-data2+(alpha.*data1))+corr;
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

c = -real(sol+x*delta_sol);
ceq = [];
end

function out = to_optimize_all(img, data1, data2,delta_sol,sol, x,normalize,with_v2, alpha)
    %Function f_global in the Masterarbeit. The subtleties (mixed 
    %minimization and normalization are considered)
    nb_elem = size(sol,1)/2;
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
    if with_v2(1)
            corr=with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
    else
        corr=0;
    end

    % Normalized minimization
    if normalize
        out=norm(Psi2w./(Psi1)-data2./(alpha.*data1))+corr;
    else
        out=norm(Psi2w-(Psi1)-data2+(alpha.*data1))+corr;
    end
    
    if isnan(out)
        out=x*10e6;
    end
end

function out = to_optimize_half(img, data1, data2,delta_sol,sol,x,normalize,with_v2,alpha)
    %Function f_half in the Masterarbeit. The subtleties (mixed 
    %minimization and normalization are considered)
    nb_elem = size(sol,1)/2;
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

    if normalize
        out=norm(Psi2w./Psi1-data2./(alpha.*data1))+corr;
    else
        out=norm(Psi2w-Psi1-data2+(alpha.*data1))+corr;        
    end

    if isnan(out)
        out=x*10e6;
    end
end
