function OptInit = init_optimization(inv_model,data1,data2,model)


nb_elem=size(inv_model.fwd_model.elems,1); % Number of finite elements of the mesh

% Initialization of the minimization vector parameter gamma_0^(tot)
try
    sol=inv_model.inv_solve_complete_diff_GN_iter.init.*ones(2*nb_elem,1);
catch
    disp('inv_solve_complete_diff_GN_iter: default initial distribution of conductivity change: 0');
    sol=(0.5+0.01*1i).*ones(2*nb_elem,1);
end

% Initialization of the img given as a result of the resolution
inv_model.jacobian_bkgnd.value=calc_cond_change(sol,nb_elem);
img = data_mapper(calc_jacobian_bkgnd( inv_model ));
img.name= 'solved by inv_solve_complete_diff_GN_iter';
img.elem_data = calc_cond_change(sol,nb_elem);
img.fwd_model= inv_model.fwd_model;

% Parameter for the step-size line-search


% Parameter for the normalization
try
    normalize=inv_model.inv_solve_complete_diff_GN_iter.normalize;
catch
    disp('inv_solve_complete_diff_GN_iter: no normalization of the measurements');
    normalize=0;
end



% Parameter to perform the alternative minimization
try
    with_v2=inv_model.inv_solve_complete_diff_GN_iter.with_v2;
    if size(with_v2,2) == 1&& with_v2==0
        with_v2 = [false,0];
    elseif size(with_v2,2) == 1&& with_v2~=0
        with_v2 = [true,with_v2];
    end
catch
    disp('inv_solve_complete_diff_GN_iter: no additional optimization regarding the measured voltages');
    with_v2=[false,0];
end




try
    data2=data2.meas;
end
try
    data1=data1.meas;
end

%% 3) Apply weighted-difference algorithm to the surfaced voltage
% Including a weighted frequency difference
try
    weighted_diff=inv_model.inv_solve_complete_diff_GN_iter.weighted_diff;
catch
    disp('inv_solve_cole_diff_GN_iter: no weigthed difference is performed');
    weighted_diff=false;
end
if weighted_diff(1)
    alpha = (dot(data1,data2)/dot(data1,data1)); % factor calculation
else
    alpha = 1;
end

opt.Display = 'iter';
opt.Maxiter = 30;
switch model
    case 'original' %only optimize conductivites of high frequency
        [step_half] = fminbnd(@(x) to_optimize_half(img,data1,data2,sol,zeros(size(sol,1),1), x,normalize,with_v2, alpha,nb_elem), ...
            1e-5, 10, opt);
        OptInit = [ones(nb_elem,1);repmat(step_half,nb_elem,1)].*sol;
    case 'absolute'
        % optimize both real and imag part of conductivites of low
        % optimum step size for low frequency, separeted for real and
        % imaginary part
        step_optim_all = fmincon(@(x) to_optimize_init_low(img, data1,sol, x,nb_elem), ...
            [1,1],[],[],[],[],[1e-5,1e-5], [10,10], [],opt);
        
        sol = complex(step_optim_all(1).*real(sol),step_optim_all(2).*imag(sol));
        
        % optimum step size for high frequency, separeted for real and
        % imaginary part
        step_optim_half=fmincon(@(x) to_optimize_init_high(img, data1, data2,sol, x,normalize,with_v2,alpha,nb_elem), ...
            [1,1],[],[],[],[],[1e-5,1e-5], [10,10],[], opt);
        sol_high = complex(step_optim_half(1).*real(sol(nb_elem+1:2*nb_elem,:)),step_optim_half(2).*imag(sol(nb_elem+1:2*nb_elem,:)));
        
        % remove imaginary part if measurement only consist of real part
        if ~any(imag(sol))
            sol = real(sol);
            sol_high = real(sol_high);
        end
        OptInit = [sol(1:nb_elem);sol_high];
        
    case 'onefunc'
        % optimize conductivites of low frequency and high frequency in
        % same time. this does not perform well - deprecated!!!
        step = fmincon(@(x) to_optimize_step(img, data1, data2,sol,zeros(size(sol,1),1), x,normalize,with_v2,alpha,nb_elem), ...
            [1,1],[],[],[],[],[1e-5,1e-5], [10,10], [],opt);
        OptInit = [repmat(step(1),nb_elem,1);repmat(step(2),nb_elem,1)].*sol;
end
end

function out = to_optimize_step(img, data1, data2,delta_sol,sol,x,normalize,with_v2,alpha,nb_elem)
img.elem_data = sol(1:nb_elem,:)+x(1)*delta_sol(1:nb_elem,:);
Psi1 = fwd_solve(img);
Psi1 = Psi1.meas;
img.elem_data = sol(nb_elem+1:2*nb_elem,:)+x(1)*x(2)*delta_sol(nb_elem+1:2*nb_elem,:);
Psi2 = fwd_solve(img);
Psi2 = Psi2.meas;
img.elem_data = alpha.*(sol(nb_elem+1:2*nb_elem,:)+x(1)*x(2)*delta_sol(nb_elem+1:2*nb_elem,:));
Psi2w = fwd_solve(img);
Psi2w = Psi2w.meas;
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


function out = to_optimize_init_low(img, data1,delta_sol, x,nb_elem)
%Function f_global in the Masterarbeit. The subtleties (mixed
%minimization and normalization are considered)

NewDelta = complex(x(1)*real(delta_sol(1:nb_elem,:)),x(2)*imag(delta_sol(1:nb_elem,:)));
img.elem_data=complex_truncation(NewDelta); %rajout de complex truncation pour eviter les conductivités negatives
Psi1=fwd_solve(img);
Psi1=Psi1.meas;

out=(norm((Psi1-data1)));

if isnan(out)
    out=x*10e6;
end
end

function out = to_optimize_init_high(img, data1, data2,delta_sol,x,normalize,with_v2,alpha,nb_elem)
%Function f_half in the Masterarbeit. The subtleties (mixed
%minimization and normalization are considered)

img.elem_data=complex_truncation(delta_sol(1:nb_elem,:));
Psi1=fwd_solve(img);
Psi1=Psi1.meas;
NewDelta = complex(x(1)*real(delta_sol(nb_elem+1:2*nb_elem,:)),x(2)*imag(delta_sol(nb_elem+1:2*nb_elem,:)));
img.elem_data=complex_truncation(NewDelta);
Psi2=fwd_solve(img);
Psi2=Psi2.meas;
img.elem_data=alpha.*complex_truncation(NewDelta);
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

function out = to_optimize_half(img, data1, data2,delta_sol,sol,x,normalize,with_v2,alpha,nb_elem)
%Function f_half in the Masterarbeit. The subtleties (mixed
%minimization and normalization are considered)

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
function cond_change=calc_cond_change(sol,nb_elem)
% Function giving the conductivity change between gamma^(high) and gamma^(low)
cond_change=sol(nb_elem+1:2*nb_elem,:)-sol(1:nb_elem,:);
end
