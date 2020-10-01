function imdl_solve = set_weighted_invprob_properties(fmdl2,init,hp_param, weighted_diff, beta,mode)
% set_weighted_invprob_properties: predefines for fdEIT reconstructions
%
% This function defines typical parameters used for the fdEIT
% reconstruction algorithm. For other specific configurations this function
% should not be used and
% imdl_solve.inv_solve_complete_diff_GN_iter. - can be set independently

% set fmdl2-Model for reconstruction
imdl_solve.fwd_model = fmdl2;
imdl_solve.type = 'inv_model';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolution (GN_complete_iterativ)                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imdl_solve.solve=@inv_solve_complete_diff_GN_iter; % Resolution using the algorithm developed for the Masterarbeit of Sébastien Dambrun

% Parameters of the resolution
imdl_solve.inv_solve_complete_diff_GN_iter.calc_step_size=true; % By default, the step size is computed, but the parameter can be changed
imdl_solve.inv_solve_complete_diff_GN_iter.nb_iter=1; % This parameter is crucil to the resolution
imdl_solve.inv_solve_complete_diff_GN_iter.init=init; % This parameter can be a scalar value or vector of length 2*number of elements of the mesh in img_solve (this in order to start with a non homogeneous background)
imdl_solve.inv_solve_complete_diff_GN_iter.normalize=false; % This to choose whether the normalized simulation must be performed or not
imdl_solve.inv_solve_complete_diff_GN_iter.display_iterations=false; % This to choose whether the voltages must be displayed at each iteration or not
imdl_solve.inv_solve_complete_diff_GN_iter.pseudo_inv=[false,1e-3]; % This to choose if the truncated Moore-Penrose pseudo-inverse must be used or not
imdl_solve.inv_solve_complete_diff_GN_iter.with_v2=beta; % This to choose if the alternative minimization must be used or not
imdl_solve.inv_solve_complete_diff_GN_iter.weighted_diff=weighted_diff; % This to choose if the weighted difference must be used or not
imdl_solve.inv_solve_complete_diff_GN_iter.bounds=[1e-6 2]; % This to give the boundaries for the line-search.
imdl_solve.inv_solve_complete_diff_GN_iter.fminbnd.Display= 'notify';
imdl_solve.inv_solve_complete_diff_GN_iter.fminbnd.MaxIter= 30;
imdl_solve.inv_solve_complete_diff_GN_iter.initial_adaptation=[false,0]; % True if the modified initial adaptation is to be performed, the integer gives the type of regularization : 0 using the mean value, 1 using the peak value.
try
    imdl_solve.inv_solve_complete_diff_GN_iter.mode = mode; % if 1, objective func= norm(dPsi-dv); if 2, objective func= norm(Psi1+Psi2w-alpha*data1-data2) 
catch 
    imdl_solve.inv_solve_complete_diff_GN_iter.mode =1;
end

imdl_solve.RtR_prior=@prior_laplace; % Choice of the prior matrix
imdl_solve.prior_laplace_noser.laplace_weight= 2; % Parameters in the case where prior_laplace_noser is chosen
imdl_solve.prior_laplace_noser.sec_order=false; % Parameters in the case where prior_laplace_noser is chosen

%    / \
%   / | \
%  /  |  \
% /   |   \
%/____o____\
% The choice of the hyperparameter is critical for the reconstruction.
% Choose wisely
imdl_solve.hyperparameter.value=hp_param; 

