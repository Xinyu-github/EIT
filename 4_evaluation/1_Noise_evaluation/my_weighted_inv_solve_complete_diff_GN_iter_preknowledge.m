function img= my_weighted_inv_solve_complete_diff_GN_iter_preknowledge( inv_model, data1, data2)
% INV_SOLVE_COMPLETE_DIFF_GN_ITER inverse solver
% img= inv_solve_complete_diff_GN_iter( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at f1
% data2      => differential data at f2
%
% both data1 and data2 may be matrices  (MxT) each of
%  M measurements at T times 
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix
%
% By default, the correct scaling of the solution that best fits the data
% is calculated at each step. 
% To provide a pre-calculated scaling (which is not recommended), specify
%     inv_model.inv_solve_complete_diff_GN_iter.calc_step_size = 0;
%     inv_model.inv_solve_complete_diff_GN_iter.step_size = [0.8,0.5]; two
%     step sizes must be specified: one is step_size, the other is
%     step_half.
% The search for correct step_size is performed using FMINBND. The default
% search interval is [1e-5 1e1]. You can modify it by specifying:
%     inv_model.inv_solve_diff_GN_one_step.bounds = [10 200];
% Additional options for FMINBD can be passed as:
%     inv_model.inv_solve_complete_diff_GN_iter.fminbnd.MaxIter = 10;
% The optimal step_size is returned in img.info.step_size.
%
% The initial conductivity distribution is 0.5+0.01i by default. It can be
% specified by the user, either as a complex number if the conductivity
% must be homogeneous, either using a vector of length 2xN (where N is the
% number of elements in the forward model) giving the initial conductivity
% of each element. This is done by specifying:
%     inv_model.inv_solve_complete_diff_GN_iter.init=desired initialization (number or vector);
%
% The default number of iterations is 1. It can be changed in the following
% way:
%     inv_model.inv_solve_complete_diff_GN_iter.nb_iter=2;
%
%
% If the alternative/mixed minimization is to be performed, the one has to
% specify it and to give the balancing parameter in the following way:
%     inv_model.inv_solve_complete_diff_GN_iter.with_v2=[true,5];
% By default, the alternative minimization is not performed.
%
% If the normalization of the input voltages must be done, it must be
% specified by the user in the following way:
%     inv_model.inv_solve_complete_diff_GN_iter.normalize=true;
% By default, the normalization is not performed.
%
% There exist the possibility to display the comparision of the input
% voltages and the voltages simulated at each iteration for the resolution.
% This is done if the user writes:
%     inv_model.inv_solve_complete_diff_GN_iter.display_iterations=true;
% By default, the iterations are not displayed.
%
% By default, a half step-size search is performed at the beginning of the
% resolution. If the modified version of this initial adaptation must be
% performed, it must be specified. There exists two types of such an
% adaptation, the one (designed by the number 0) uses the average value of 
% the voltage vectors, the other (designed by the number 1) uses the
% maximal value. The way this parameter must be specified is:
%     inv_model.inv_solve_complete_diff_GN_iter.initial_adaptation=[true,0];
%
% For stability reasons, it may be possible that a truncated Moore-Penrose
% pseudo-inverse of the reconstruction matrix must be used instead of its 
% inverse. This can be decidecd by the user with:
%    inv_model.inv_solve_complete_diff_GN_iter.pseudo_inv=[true;1e-5];
% where the second parameter gives the threshold for the truncation of the
% eigenvalues (see function trunc_pseudo_inv at the end of this file)
%
% See also INV_SOLVE, CALC_SOLUTION_ERROR, FMINBND
clear ss ii vals;
global ii vals ss;
ii = 1; vals = 1; ss(ii) = vals; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Initialization of all the reconstruction parameters           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_elem=size(inv_model.fwd_model.elems,1); % Number of finite elements of the mesh

% Initialization of the minimization vector parameter gamma_0^(tot)
try
    sol=(inv_model.my_weighted_inv_solve_complete_diff_GN_iter.init).*ones(2*nb_elem,1);
    ii = ii+1; vals = 2; ss(ii) = vals; 
catch
    disp('my_weighted_inv_solve_complete_diff_GN_iter: default initial distribution of conductivity change: 0');
    sol=(0.5+0.01*1i).*ones(2*nb_elem,1);
    ii = ii+1; vals = 3; ss(ii) = vals; 
end

% Initialization of the img given as a result of the resolution
gamma_diff = calc_cond_change(sol,nb_elem); 
inv_model.jacobian_bkgnd.value=gamma_diff;
img = data_mapper(calc_jacobian_bkgnd( inv_model ));
img.name= 'solved by my_weighted_inv_solve_complete_diff_GN_iter';
img.elem_data = calc_cond_change(sol,nb_elem);
img.fwd_model= inv_model.fwd_model;

% Parameter for the step-size line-search
try 
 opt = inv_model.my_weighted_inv_solve_complete_diff_GN_iter.fminbnd;
     ii = ii+1; vals = 4; ss(ii) = vals; 
catch
 opt.Display = 'iter';
 ii = ii+1; vals = 5; ss(ii) = vals;
end

% Parameter for the normalization
try
    normalize=inv_model.my_weighted_inv_solve_complete_diff_GN_iter.normalize;
    ii = ii+1; vals = 6; ss(ii) = vals;
catch
    disp('my_weighted_inv_solve_complete_diff_GN_iter: no normalization of the measurements');
    normalize=0;
    ii = ii+1; vals = 7; ss(ii) = vals;
end

% Parameter to display the voltages at each iterations
try
    display_iterations=inv_model.my_weighted_inv_solve_complete_diff_GN_iter.display_iterations;
    ii = ii+1; vals = 8; ss(ii) = vals;
catch
    disp('inv_solve_complete_diff_GN_iter: no display of the iteration results');
    display_iterations=false;
    ii = ii+1; vals = 9; ss(ii) = vals;
end

% Parameter to perform the alternative minimization
try
    with_v2=inv_model.my_weighted_inv_solve_complete_diff_GN_iter.with_v2;
    ii = ii+1; vals = 10; ss(ii) = vals;
catch
    disp('my_weighted_inv_solve_complete_diff_GN_iter: no additional optimization regarding the measured voltages');
    with_v2=[false,1];
    ii = ii+1; vals = 11; ss(ii) = vals;
end

% Parameter to perform the modified initial adaptation
try
    init_adapt=inv_model.my_weighted_inv_solve_complete_diff_GN_iter.initial_adaptation;
    ii = ii+1; vals = 12; ss(ii) = vals;
catch
    disp('my_weighted_inv_solve_complete_diff_GN_iter: no adaptation of the initial value');
    init_adapt=[false,0];
    ii = ii+1; vals = 13; ss(ii) = vals;
end

% Parameter giving the total number of iterations
try
    nb_iter=inv_model.my_weighted_inv_solve_complete_diff_GN_iter.nb_iter;
    ii = ii+1; vals = 14; ss(ii) = vals;
catch
    disp('my_weighted_inv_solve_complete_diff_GN_iter: default number of iterations: 1');
    nb_iter=1;
    ii = ii+1; vals = 15; ss(ii) = vals;
end

% Parameter to use the truncated Moore-Penrose pseudo-inverse instead of
% the inverse of the reconstruction matrix. See the function
% trunc_pseudo_inv at the end of this file for more infos
try
    pinv=inv_model.my_weighted_inv_solve_complete_diff_GN_iter.pseudo_inv;
    ii = ii+1; vals = 16; ss(ii) = vals;
catch
    disp('my_weighted_inv_solve_cole_diff_GN_iter_bis: calculation of the inverse');
    pinv=[false,1];
    ii = ii+1; vals = 17; ss(ii) = vals;
end

% Initialization of the hyperparameter
hp  = calc_hyperparameter(inv_model); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Start of the resolution                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial adapation of the solution. Either half line search or
% multiplication by a chosen constant
if init_adapt(1)    
    sol=initial_adaptation(sol, data1, data2, inv_model, init_adapt(2));
    ii = ii+1; vals = 18; ss(ii) = vals;
else
alpha = inv_model.my_weighted_inv_solve_complete_diff_GN_iter.weighted_diff;
step_optim=fminbnd(@(x) to_optimize_half(img, data1, data2,sol,zeros(size(sol,1),1), x,normalize,with_v2,alpha), ...
                           1e-5, 10, opt);
sol=[ones(nb_elem,1);step_optim*ones(nb_elem,1)].*sol;
ii = ii+1; vals = 19; ss(ii) = vals;
end

% Normalization of the differential input data if it was so decided %ss(63)
if normalize
    try
        dv = calc_difference_data(data1,data2, inv_model.fwd_model)./data1;
        ii = ii+1; vals = 20; ss(ii) = vals;
    catch
        dv = calc_difference_data(data1,data2, inv_model.fwd_model)./data1.meas;
        ii = ii+1; vals = 21; ss(ii) = vals;
    end
else
       % delta U
        dv = calc_difference_data(data1,data2, inv_model.fwd_model);
        ii = ii+1; vals = 22; ss(ii) = vals;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Start of the iterative minimization                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ss(88)-ss()
alpha = inv_model.my_weighted_inv_solve_complete_diff_GN_iter.weighted_diff;
for i=1:nb_iter
    ii = ii+1; vals = 23; ss(ii) = vals;
%% calulate RM(Hessian),RtR(prior matrix) + arrange J(jacobi matrix) using "get_RM"    
    % The reconstruction matrix is either calculated or retrieved from the
    % cache. The function for the calculation of the reconstruction matrix
    % is get_RM, it is written in this file
    [RM,J,L]= eidors_cache(@get_RM, {inv_model,hp,sol,nb_elem,normalize,with_v2,alpha});
    img_fwd_model=mk_image(inv_model);
     save('L','L');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate delta A = A(cond2) - alpha A(cond1) and manage U1, U2
    % Simulation of the voltage measurement with the current value of the
    % conductivity vector/minimization variable gamma_k^(tot)
    img_fwd_model.elem_data=sol(1:nb_elem,:);
    Psi1=fwd_solve(img_fwd_model);
    Psi1=Psi1.meas;

    img_fwd_model.elem_data=sol((nb_elem+1):2*nb_elem,:);
    Psi2=fwd_solve(img_fwd_model);
    Psi2=Psi2.meas;

    % If necessary, the input data are converted from the EIDORS-data
    % datatype to a vector.
    try
        data2=data2.meas;
        ii = ii+1; vals = 24; ss(ii) = vals;
    end
    
    try
        data1=data1.meas;
        ii = ii+1; vals = 24; ss(ii) = vals;
    end
    
    % If the normalization is to be performed, then the simulated data are
    % also normalized.
    switch normalize
        case 0
            ii = ii+1; vals = 25; ss(ii) = vals;
            if norm(Psi1) < norm(Psi2)
                
                
                dPsi = Psi2 - alpha.*Psi1;  
                
                
            else
%                 dPsi = alpha.*Psi2 - Psi1;
                  dPsi = Psi2 - alpha.*Psi1;
            end
        case 1
            dPsi = (Psi2 - Psi1)./Psi1;
            ii = ii+1; vals = 26; ss(ii) = vals;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate inv of RM or H^{-1} 
    % If the truncated pseudo-inverse must be used then it is calculated
    if pinv(1)
        RM=trunc_pseudo_inv(RM,pinv(2));
        ii = ii+1; vals = 27; ss(ii) = vals;
    else
        RM=inv(RM);
        ii = ii+1; vals = 28; ss(ii) = vals;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate delta conductivity (and grad)
    % Calculation of the update step. Two possibilities, depending on the
    % problem formulation (with oder without mixed minimization)
    if with_v2(1)
        %delta_sol = -H^(-1)* grad 
        delta_sol = -RM     *( J'*[ dPsi-dv; 
                                    with_v2(2)*(Psi1-data1); 
                                    with_v2(2)*(Psi2-data2)  ] + hp*L*sol );
        ii = ii+1; vals = 29; ss(ii) = vals;
    else
        delta_sol = -RM*( J'*(dPsi-dv) + hp*L*sol );
        ii = ii+1; vals = 30; ss(ii) = vals;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % If the voltages must be displayed, it is done here
    if display_iterations
        plot_simulation_data([2*(i-1)+1,2*i],sol, delta_sol, Psi1, Psi2, data1, data2,i-1,normalize==1);
        ii = ii+1; vals = 31; ss(ii) = vals;
    end
    
%% Calculate optimization step size (α_global, α_half) from 2-line search using scale_to_fit_data
    % The two line-search are here performed
    [step_all,step_half] = scale_to_fit_data( img_fwd_model,inv_model, data1, data2,delta_sol,sol,with_v2);

%% Update conductivity solution conductivity_k+1
    % The function complex_truncation ensures that there are no
    % conductivities with a negative real part (which is not physically possible)
    sol=complex_truncation(sol+[step_all*ones(nb_elem,1);step_all*step_half*ones(nb_elem,1)].*delta_sol);

end

% If the voltages must be displayed at each iteration, it is done after the
% en dof the loop
if display_iterations
    ii = ii+1; vals = 32; ss(ii) = vals;
    img_fwd_model.elem_data=sol(1:nb_elem,:);
    Psi1=fwd_solve(img_fwd_model);
    Psi1=Psi1.meas;


    img_fwd_model.elem_data=sol((nb_elem+1):2*nb_elem,:);
    Psi2=fwd_solve(img_fwd_model);
    Psi2=Psi2.meas;
    
    plot_simulation_data([2*nb_iter+1,2*nb_iter+2],sol, delta_sol, Psi1, Psi2, data1, data2,i,normalize==1)
end

%% End of the reconstruction, the resulting image is updated with the correct values.
img.elem_data=calc_cond_change(sol,nb_elem);
img.cond1=sol(1:nb_elem);
img.cond2=sol(nb_elem+1:end);

imgr = load('img_diff'); imgr = imgr.img_diff;


elem = img.elem_data-min(img.elem_data);
elem = elem/max(elem);

elem = elem.*(max(abs(imgr.elem_data)));
elem = elem+min(imgr.elem_data);

img.elem_data = elem;

save('sol','sol'); save('Psi1','Psi1'); save('Psi2','Psi2'); save('nb_elem','nb_elem');
save('to_check_inv_w','ss'); 
end
   
        

       