function img= inv_solve_complete_diff_GN_iter_mixterm( inv_model, data1, data2)
% INV_SOLVE_COMPLETE_DIFF_GN_ITER inverse solver
% img= inv_solve_complete_diff_GN_iter( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at f1
% data2      => differential data at f2
%
% both data1 and data2 may be matrices (MxT) each of
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Initialization of all the reconstruction parameters           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
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
try 
 opt = inv_model.inv_solve_complete_diff_GN_iter.fminbnd;
catch
 opt.Display = 'iter';
end
mix_=inv_model.inv_solve_complete_diff_GN_iter.mix_scale;
% Parameter for the normalization
try
    normalize=inv_model.inv_solve_complete_diff_GN_iter.normalize;
catch
    disp('inv_solve_complete_diff_GN_iter: no normalization of the measurements');
    normalize=0;
end

% Parameter to display the voltages at each iterations
try
    display_iterations=inv_model.inv_solve_complete_diff_GN_iter.display_iterations;
catch
    disp('inv_solve_complete_diff_GN_iter: no display of the iteration results');
    display_iterations=false;
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
    with_v2=[false,1];
end

% Parameter to perform the modified initial adaptation
try
    init_adapt=inv_model.inv_solve_complete_diff_GN_iter.initial_adaptation;
catch
    disp('inv_solve_complete_diff_GN_iter: no adaptation of the initial value');
    init_adapt=[false,0];
end

% Parameter giving the total number of iterations
try
    nb_iter=inv_model.inv_solve_complete_diff_GN_iter.nb_iter;
catch
    disp('inv_solve_complete_diff_GN_iter: default number of iterations: 1');
    nb_iter=1;
end

% Parameter to use the truncated Moore-Penrose pseudo-inverse instead of
% the inverse of the reconstruction matrix. See the function
% trunc_pseudo_inv at the end of this file for more infos
try
    pinv=inv_model.inv_solve_complete_diff_GN_iter.pseudo_inv;
catch
    disp('inv_solve_cole_diff_GN_iter_bis: calculation of the inverse');
    pinv=[false,1];
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

% Initialization of the hyperparameter
hp  = calc_hyperparameter(inv_model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Start of the resolution                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial adapation of the solution. Either half line search or
% multiplication by a chosen constant
if init_adapt(1)    
    sol=initial_adaptation(sol, data1, data2, inv_model, init_adapt(2), alpha);
else
% step_optim=fminbnd(@(x) to_optimize_half(img, data1, data2,sol,zeros(size(sol,1),1), x,normalize,with_v2,alpha), ...
%                            1e-5, 10, opt);
step_optim=1;
sol=[ones(nb_elem,1);step_optim*ones(nb_elem,1)].*sol;
end

% Normalization of the differential input data if it was so decided
if normalize
    try
        dv = calc_difference_data(alpha.*data1,data2, inv_model.fwd_model)./(alpha.*data1);
    end
else
        dv = calc_difference_data(alpha.*data1,data2, inv_model.fwd_model);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Start of the iterative minimization                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Laufzeit Initialization:')
toc
for i=1:nb_iter
    tic
    % The reconstruction matrix is either calculated or retrieved from the
    % cache. The function for the calculation of the reconstruction matrix
    % is get_RM, it is written in this file
    [RM,J,LtL,L]= eidors_cache(@get_fdRM, {inv_model,hp,sol,nb_elem,normalize,with_v2,alpha,mix_});
    img_fwd_model=mk_image(inv_model);

    % Simulation of the voltage measurement with the current value of the
    % conductivity vector/minimization variable gamma_k^(tot)
    img_fwd_model.elem_data=sol(1:nb_elem,:);
    Psi1=fwd_solve(img_fwd_model);
    Psi1=Psi1.meas;

    img_fwd_model.elem_data=sol((nb_elem+1):2*nb_elem,:);
    Psi2=fwd_solve(img_fwd_model);
    Psi2=Psi2.meas;
    
    img_fwd_model.elem_data=alpha.*sol((nb_elem+1):2*nb_elem,:);
    Psi2w=fwd_solve(img_fwd_model);
    Psi2w=Psi2w.meas;

    % If necessary, the input data are converted from the EIDORS-data
    % datatype to a vector.
    try
        data2=data2.meas;
    end
    
    try
        data1=data1.meas;
    end
    
    % If the normalization is to be performed, then the simulated data are
    % also normalized.
    switch normalize
        case 0
            dPsi=Psi2w-Psi1;
        case 1
            dPsi=(Psi2w-Psi1)./Psi2w;
    end

    % If the truncated pseudo-inverse must be used then it is calculated
%     if pinv(1)
%         RM=trunc_pseudo_inv(RM,pinv(2));
%     else
%         RM=inv(RM);
%     end

    % Calculation of the update step. Two possibilities, depending on the
    % problem formulation (with oder without mixed minimization)
    if with_v2(1)
        delta_sol=-RM\(J'*[dPsi-dv;with_v2(2)*(Psi1-data1);with_v2(2)*(Psi2-data2);mix_*(Psi1.*Psi2-data1.*data2)]+hp*L'*sol);
    else
        delta_sol=-RM\(J'*[dPsi-dv;mix_*(Psi1.*Psi2-data1.*data2)]+hp*L'*sol);
    end
    
    % If the voltages must be displayed, it is done here
    if display_iterations
        plot_simulation_data([2*(i-1)+1,2*i],sol, delta_sol, Psi1, Psi2, data1, data2,i-1,normalize==1);
    end

    % The two line-search are here performed
    [step_all,step_half] = scale_to_fit_data( img_fwd_model,inv_model, data1, data2,delta_sol,sol,with_v2,alpha,mix_);
    % The function complex_truncation ensures that there are no
    % conductivities with a negative real part (which is not physically possible)
    sol=complex_truncation(sol+[step_all*ones(nb_elem,1);step_all*step_half*ones(nb_elem,1)].*delta_sol);
    disp(strcat('Laufzeit Iteration ',num2str(i),':'))
    toc
end

% If the voltages must be displayed at each iteration, it is done after the
% en dof the loop
if display_iterations
    img_fwd_model.elem_data=sol(1:nb_elem,:);
    Psi1=fwd_solve(img_fwd_model);
    Psi1=Psi1.meas;


    img_fwd_model.elem_data=sol((nb_elem+1):2*nb_elem,:);
    Psi2=fwd_solve(img_fwd_model);
    Psi2=Psi2.meas;
    
    plot_simulation_data([2*nb_iter+1,2*nb_iter+2],sol, delta_sol, Psi1, Psi2, data1, data2,i,normalize==1)
end

% End of the reconstruction, the resulting image is updated with the
% correct values.
img.elem_data=calc_cond_change(sol,nb_elem);
img.cond1=sol(1:nb_elem);
img.cond2=sol(nb_elem+1:end);
    img_fwd_model.elem_data=sol(1:nb_elem,:);
    Psi1=fwd_solve(img_fwd_model);
img.Psi1=Psi1.meas;
    img_fwd_model.elem_data=sol((nb_elem+1):2*nb_elem,:);
    Psi2=fwd_solve(img_fwd_model);
img.Psi2=Psi2.meas;
    img_fwd_model.elem_data=alpha*sol((nb_elem+1):2*nb_elem,:);
    Psi2w=fwd_solve(img_fwd_model);
img.Psi2w=Psi2w.meas;
img.data1=data1;
img.data2=data2;
img.inv_model = inv_model;
img.L=L;
img.alpha=alpha;
img.hp=hp;
img.beta = with_v2(2); 
% img.step_optim = step_optim;
img.step_size = [step_all,step_half];
end



function [RM,J,RtR,R] = get_fdRM( inv_model,hp,bkgnd,nb_elem,normalize,with_v2,alpha,mix_)
% This function is used to compute the reconstruction matrix RM, the
% Jacobian of the objective function J and the prior matrix RtR 
    inv_model.jacobian_bkgnd.value=bkgnd(1:nb_elem);
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
% If the normalized version is to be used, the voltages measurements must
% be used in the jacobian (cf schriftliche Fassung der Masterarbeit Seite 34)
    
        Psi1=fwd_solve(img_bkgnd);
        Psi1=Psi1.meas;
   
   % Calculation of the sensitivity matrix for the low-frequency
   % conductivity
   J1 = calc_jacobian( img_bkgnd);  
   
   % Same thing for the high-frequency conductivity
   inv_model.jacobian_bkgnd.value=alpha.*bkgnd((nb_elem+1):2*nb_elem);
   img_bkgnd= calc_jacobian_bkgnd( inv_model );
 
   
   J2 = calc_jacobian( img_bkgnd);
   inv_model.jacobian_bkgnd.value=bkgnd((nb_elem+1):2*nb_elem);
   img_bkgnd= calc_jacobian_bkgnd( inv_model );
   
        Psi2=fwd_solve(img_bkgnd);
        Psi2=Psi2.meas;
 
   J3 = calc_jacobian( img_bkgnd);
 
   % The Jacobian of the total problem is defined. It is different if we
   % perform the normalized minimization (cf Masterarbeit Seite 34)
    if normalize==1
        J=[-Psi2.*J1./(Psi1.^2),J2./Psi1];
    else
        J=[-J1,J2];
    end
    % The Jacobian must also be adapted to the mixed minimization (cf Masterarbeit Seite 31)
    if with_v2(1)
        J=[J;with_v2(2)*J1,zeros(size(J1));zeros(size(J1)),with_v2(2)*J2];
    end
    J = [J;mix_.*[J1.*Psi2,J3.*Psi1]];
    % Calculation of the desired prior-matrix with the functiom calc_RtR_prior 
   inv_model.jacobian_bkgnd.value=real(bkgnd(1:nb_elem));
%    RtR = calc_RtR_prior( inv_model );
    R = calc_R_prior(inv_model);
%    RtR=RtR'*RtR;
    
     RtR = R;
      R =blkdiag(R,R);
   %The prior matrix is adapted to the formulation of the problem with
   %gamma_0^(tot) (cf Masterarbeit Seite 27)
   RtR=blkdiag(RtR,RtR);
   J_sym=J'*J;
    %JtJ results in mixed gradients of Jhigh and Jlow -> set to zero
    J_sym(1:nb_elem,(nb_elem+1):2*nb_elem) = 0;
    J_sym((nb_elem+1):2*nb_elem,1:nb_elem) = 0;
   
   % Definition of the reconstruction matrix
   RM=J_sym + hp*RtR;
end

   
   
function [step_all,step_half] = scale_to_fit_data(img, inv_model,data1, data2,delta_sol,sol, with_v2, alpha,mix_)
   % Function to find the step size to multiply sol by to best fit data
   step_all = 1;
   step_half = 1;
   % By default, the step size is computed.
   do_step   = true;
   try do_step = inv_model.inv_solve_complete_diff_GN_iter.calc_step_size; 
   end
   
   if do_step
       %if the step-size must be calculated
      eidors_msg('inv_solve_complete_diff_GN_iter: Calculating optimal step size to fit data',2);
      % options for fminbnd
      try 
         opt = inv_model.inv_solve_complete_diff_GN_iter.fminbnd;
      catch
         opt.Display = 'iter';
         opt.MaxIter = 30 ;
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
        
        %In order not to have negative conductivities after the update, we
        %must enforce a condition on the step-size : it must not be
        %negative (which is obvious considering the Gauss-Newton algorithm)
        %and it must not be larger than the value x_range2.
        
        aux=-real(sol)./real(delta_sol);
        x_range2=min(aux(aux>range(1)));
        
        if isempty(x_range2) 
            x_range2=range(2);
        elseif x_range2>range(2)
            x_range2=range(2);
        end
        
        %First step-size line-search corresponding to f_global in the
        %Masterabeit (Seite 27)
       clearvars x;
      [step_all,fval1] = fminbnd(@(x) to_optimize_all(img,data1,data2,delta_sol,sol, x,normalize,with_v2, alpha,mix_), ...
                            range(1), x_range2, opt);
                        
       nb_elem=size(sol)/2;
       
       % Same precaution as before in order not to have conductivities with
       % a negative real part.
       aux=-real(sol(nb_elem+1:end))./(step_all*real(delta_sol(nb_elem+1:end)));
       x_range2=min(aux(aux>range(1)));
       
        if isempty(x_range2) 
            x_range2=range(2);
        elseif x_range2>range(2)
            x_range2=range(2);
        end
           
         %Second step-size line-search corresponding to f_half in the
        %Masterabeit (Seite 28)
        clearvars x;
      [step_half,fval2] = fminbnd(@(x) to_optimize_half(img,data1,data2,step_all*delta_sol,sol, x,normalize,with_v2, alpha,mix_), ...
                              range(1), x_range2, opt);
      % If the second step-size line-search didn't improve anything, then
      % ignore it
      if fval2>fval1
          step_half=1;
      end
                          
   else
      % if not calculating, check if step_size provided
      try
         [step_all,step_half] = inv_model.inv_solve_complete_diff_GN_one_step.step_size;
      end
   end
   nb_elem=size(sol,1)/2;

   % complex_truncation enforces a last time the condition that no
   % conductivity has a negative real part
   img.elem_data = complex_truncation(calc_cond_change(sol+[step_all*ones(nb_elem,1);step_all*step_half*ones(nb_elem,1)].*delta_sol,nb_elem));
   try 
       img.info.step_size=[img.info.step_size;step_all,step_half];
   catch
       img.info.step_size = [step_all,step_half];
   end
end

function out = to_optimize_all(img, data1, data2,delta_sol,sol, x,normalize,with_v2, alpha,mix_)
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
    if with_v2(1)
            corr=with_v2(2)*(norm((Psi1-data1))+norm((Psi2-data2)));
    else
        corr=0;
    end

    % Normalized minimization
    if normalize
        out=norm(Psi2w./(Psi1)-data2./(alpha.*data1))+corr;
    else
        out=norm(Psi2w-(Psi1)-data2+(alpha.*data1))+corr+mix_*norm(Psi1.*Psi2-data1.*data2);
    end
    
    if isnan(out)
        out=x*10e6;
    end
end
    
    function out = to_optimize_half(img, data1, data2,delta_sol,sol,x,normalize,with_v2,alpha,mix_)
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

   	if normalize
        out=norm(Psi2w./Psi1-data2./(alpha.*data1))+corr;
    else
        out=norm(Psi2w-Psi1-data2+(alpha.*data1))+corr+mix_*norm(Psi1.*Psi2-data1.*data2);        
    end
    
    if isnan(out)
        out=x*10e6;
    end
    end

   
function cond_change=calc_cond_change(sol,nb_elem)
% Function giving the conductivity change between gamma^(high) and gamma^(low)
cond_change=sol(nb_elem+1:2*nb_elem,:)-sol(1:nb_elem,:);
end

function new_sol=initial_adaptation(sol, data1, data2, inv_model, adapt_type, alpha)
% Function to perform the modified initial adaptation (Seite 32 der Masterarbeit)
    try
        data2=data2.meas;
    end   
    try
        data1=data1.meas;
    end
    
    if adapt_type>1 || adapt_type<0
        disp('The initial adaptation type must be either 0 or 1. No initial adaptation was performed');
        new_sol=sol;
        return
    end
    
    nb_elem=length(sol)/2;
    img_adapt=mk_image(inv_model);
    img_adapt.elem_data=sol(1:nb_elem);
    Psi1=fwd_solve(img_adapt);
    Psi1=(Psi1.meas);
    
    img_adapt.elem_data=alpha.*sol(nb_elem+1:end);
    Psi2=fwd_solve(img_adapt);
    Psi2=Psi2.meas;
    
    switch adapt_type
        case 0
            a1=mean(data1)/mean(Psi1);
            a2=mean(data2)/mean(Psi2);


        case 1

            a1=max(data1)/max(Psi1);
            a2=max(data2)/max(Psi2);
        
        otherwise
            error('inv_solve_complete_diff_GN_iter: adaptation type not valid');
               
    end
    
    new_sol=sol(1:nb_elem)/a1;
    new_sol=[new_sol;sol(nb_elem+1:end)/a2];
    
    
end

function A=trunc_pseudo_inv(M,tol)
% Function giving the truncated Moore-Penrose pseudo-inverse of a matrix.
% Only the biggest "eigenvalues" are conserved: all the eigenvalues bigger
% than tol*(max(eigenvalues)). This so that the smaller eigenvalues do not
% introduce instabilities in the pseudo-inverse.

    [U,S,V]=svd(M,'econ');
    singular_values=diag(S);
    singular_values=singular_values.*(singular_values>tol*singular_values(1));
    for idx=1:min(size(S))
        if singular_values(idx)~=0
        S(idx,idx)=1/singular_values(idx);
        end
    end

    A=V*S'*U';
end
   
        

       
