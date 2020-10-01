function [RM,J,RtR] = get_fd_weightedRM( inv_model,bkgnd)
% This function is used to compute the reconstruction matrix RM, the
% Jacobian of the objective function J and the prior matrix RtR 
    try
        with_v2=inv_model.inv_solve_complete_diff_GN_iter.with_v2;
    catch
        disp('inv_solve_complete_diff_GN_iter: no additional optimization regarding the measured voltages');
        with_v2=[false,0];
    end
  
    hp = calc_hyperparameter( inv_model );
    
    inv_model.jacobian_bkgnd.value=1.2*bkgnd;
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
% If the normalized version is to be used, the voltages measurements must
% be used in the jacobian (cf schriftliche Fassung der Masterarbeit Seite 34)
   
   % Calculation of the sensitivity matrix for the low-frequency
   % conductivity
   J1 = calc_jacobian( img_bkgnd);  
   % Same thing for the high-frequency conductivity
   inv_model.jacobian_bkgnd.value=1.5*bkgnd;
   img_bkgnd= calc_jacobian_bkgnd( inv_model );

   J2 = calc_jacobian( img_bkgnd);
 
   % The Jacobian of the total problem is defined. It is different if we
   % perform the normalized minimization (cf Masterarbeit Seite 34)

   J=[-J1,J2];
    
    % The Jacobian must also be adapted to the mixed minimization (cf Masterarbeit Seite 31)
    if with_v2(1)
        J=[J;with_v2(2)*J1,zeros(size(J1));zeros(size(J1)),with_v2(2)*J2];
        J_sym=J'*J;
        W = ones(size(J_sym)).*((1/(1+with_v2(2)^2))*blkdiag(ones(size(J_sym)/2),ones(size(J_sym)/2)));
    else
        J_sym=J'*J;
        W = ones(size(J_sym));
    end
    
    % Calculation of the desired prior-matrix with the functiom calc_RtR_prior 
   inv_model.jacobian_bkgnd.value=bkgnd;
   RtR = calc_RtR_prior( inv_model );

%    RtR=RtR'*RtR;
    
   %The prior matrix is adapted to the formulation of the problem with
   %gamma_0^(tot) (cf Masterarbeit Seite 27)
   RtR=blkdiag(RtR,RtR);
   
    

   % Definition of the reconstruction matrix
   RM=W.*J_sym + hp*RtR;
end