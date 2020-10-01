function imdl_solve = mk_hpselection_imdl(fmdl,imdl_solve,hp,beta,method,init,weighted_diff,mode)
% set fmdl-Model for reconstruction
imdl_solve.fwd_model = fmdl;
imdl_solve.rec_model = imdl_solve;
imdl_solve.type = 'inv_model';
imdl_solve.hyperparameter.value=hp; 
if beta ~= 0
    imdl_solve.inv_solve_complete_diff_GN_iter.with_v2=[true,beta]; % This to choose if the alternative minimization must be used or not
else
    imdl_solve.inv_solve_complete_diff_GN_iter.with_v2=[false,0]; % This to choose if the alternative minimization must be used or not   
end
imdl_solve.RtR_prior=@prior_laplace; % Choice of the prior matrix
imdl_solve.prior_laplace_noser.laplace_weight= 2; % Parameters in the case where prior_laplace_noser is chosen
imdl_solve.prior_laplace_noser.sec_order=false; % Parameters in the case where prior_laplace_noser is chosen
imdl_solve.inv_solve_complete_diff_GN_iter.weighted_diff=weighted_diff; % This to choose if the weighted difference must be used or not
 imdl_solve.inv_solve_complete_diff_GN_iter.normalize=false; % This to choose whether the normalized simulation must be performed or not
imdl_solve.inv_solve_complete_diff_GN_iter.mode=mode;               
switch method
    case 'lcurve'
        imdl_solve.solve=@inv_solve_complete_diff_GN_iter; % Resolution using the algorithm developed for the Masterarbeit of Sébastien Dambrun
        imdl_solve.inv_solve_complete_diff_GN_iter.calc_step_size=true; % By default, the step size is computed, but the parameter can be changed
        imdl_solve.inv_solve_complete_diff_GN_iter.nb_iter=1; % This parameter is crucil to the resolution
        imdl_solve.inv_solve_complete_diff_GN_iter.init=init; % This parameter can be a scalar value or vector of length 2*number of elements of the mesh in img_solve (this in order to start with a non homogeneous background)
        imdl_solve.inv_solve_complete_diff_GN_iter.display_iterations=false; % This to choose whether the voltages must be displayed at each iteration or not
        imdl_solve.inv_solve_complete_diff_GN_iter.pseudo_inv=[false,1e-3]; % This to choose if the truncated Moore-Penrose pseudo-inverse must be used or not
        imdl_solve.inv_solve_complete_diff_GN_iter.step_optim =0;

        imdl_solve.inv_solve_complete_diff_GN_iter.bounds=[1e-6 2]; % This to give the boundaries for the line-search.
        imdl_solve.inv_solve_complete_diff_GN_iter.initial_adaptation=[false,0]; % True if the modified initial adaptation is to be performed, the integer gives the type of regularization : 0 using the mean value, 1 using the peak value.

        imdl_solve.RtR_prior=@prior_laplace; % Choice of the prior matrix
        imdl_solve.prior_laplace_noser.laplace_weight= 2; % Parameters in the case where prior_laplace_noser is chosen
        imdl_solve.prior_laplace_noser.sec_order=false; % Parameters in the case where prior_laplace_noser is chosen


    case 'bestres'
        nb_noise = 50;% Specify the number of noise
        imdl_solve.NoiseSNR = linspace(1e3,5e3,nb_noise);% Specify noise amplitude
    case 'snr'
        imdl_solve.hyperparameter.n_targets = 700;
        
        

end

