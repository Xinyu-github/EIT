function fmdl = my_weighted_invprob_properties(fmdl2,init_iter,hp_param,weighted_diff,beta)
% fmdl2 = load('imdl_solve_sim');
% fmdl.fwd_model = fmdl2.imdl_solve_sim.fwd_model;
% fmdl2 = load('imdl_solve_sim7');
fmdl.fwd_model = fmdl2;
fmdl.type = 'inv_model';
% fmdl.gnd_node = 1;
fmdl.solve=@my_weighted_inv_solve_complete_diff_GN_iter; % Resolution using the algorithm developed for the Masterarbeit of S�bastien Dambrun

% Parameters of the resolution
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.calc_step_size=true; % By default, the step size is computed, but the parameter can be changed
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.nb_iter=2; % This parameter is crucil to the resolution
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.init=init_iter;%0.05+1i*0.03;%0.01+1i*0.01; % This parameter can be a scalar value or vector of length 2*number of elements of the mesh in img_solve (this in order to start with a non homogeneous background)
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.normalize=false; % This to choose whether the normalized simulation must be performed or not
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.display_iterations=false; % This to choose whether the voltages must be displayed at each iteration or not
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.pseudo_inv=[false,1e-3]; % This to choose if the truncated Moore-Penrose pseudo-inverse must be used or not
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.with_v2=[true,beta]; % This to choose if the alternative minimization must be used or not
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.bounds=[1e-6 2]; % This to give the boundaries for the line-search.
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.initial_adaptation=[false,0]; % True if the modified initial adaptation is to be performed, the integer gives the type of regularization : 0 using the mean value, 1 using the peak value.
fmdl.my_weighted_inv_solve_complete_diff_GN_iter.weighted_diff=weighted_diff;

fmdl.RtR_prior=@prior_laplace; % Choice of the prior matrix
fmdl.prior_laplace_noser.laplace_weight= 2; % Parameters in the case where prior_laplace_noser is chosen
fmdl.prior_laplace_noser.sec_order=false; % Parameters in the case where prior_laplace_noser is chosen

%    / \
%   / | \
%  /  |  \
% /   |   \
%/____o____\
% The choice of the hyperparameter is critical for the reconstruction.
% Choose wisely
fmdl.hyperparameter.value=hp_param;%10^(-6);%10^(-4); 