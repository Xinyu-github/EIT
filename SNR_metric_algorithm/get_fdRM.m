function RM = get_fdRM( inv_model,bkgnd1,bkgnd2,alpha)
% This function is used to compute the reconstruction matrix RM, the
% Jacobian of the objective function J and the prior matrix RtR 
    nb_elem = size(inv_model.fwd_model.elems,1);
    hp = calc_hyperparameter(inv_model);
    try
        mode = inv_model.inv_solve_complete_diff_GN_iter.mode;
    catch
        mode = 1;
    end
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
    if mode == 2
        with_v2 = [1,1];
    end
    inv_model.jacobian_bkgnd.value=bkgnd1;
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    % If the normalized version is to be used, the voltages measurements must
    % be used in the jacobian (cf schriftliche Fassung der Masterarbeit Seite 34)
%     if normalize==1
%         Psi1=fwd_solve(img_bkgnd);
%         Psi1=Psi1.meas;
%         Psi1=repmat(Psi1,1,nb_elem);
%     end
%     % Calculation of the sensitivity matrix for the low-frequency
    % conductivity
    J1 = calc_jacobian( img_bkgnd);  

    % Same thing for the high-frequency conductivity
    inv_model.jacobian_bkgnd.value=alpha.*bkgnd2;
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
%     if normalize==1
%     Psi2=fwd_solve(img_bkgnd);
%     Psi2=Psi2.meas;
%     Psi2=repmat(Psi2,1,nb_elem);
%     end
    J2 = calc_jacobian( img_bkgnd);

    inv_model.jacobian_bkgnd.value=bkgnd2;
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
%     if normalize==1
%     Psi2=fwd_solve(img_bkgnd);
%     Psi2=Psi2.meas;
%     Psi2=repmat(Psi2,1,nb_elem);
%     end
    J3 = calc_jacobian( img_bkgnd);

    % The Jacobian of the total problem is defined. It is different if we
    % perform the normalized minimization (cf Masterarbeit Seite 34)
%     if normalize==1
%         if mode ==1
%             J=[-Psi2.*J1./(Psi1.^2),J2./Psi1;with_v2(2)*J1,zeros(size(J1));zeros(size(J1)),with_v2(2)*J2];
%         elseif mode == 2
%             J=[-Psi2.*J1./(Psi1.^2),J2./Psi1;J1,zeros(size(J1));zeros(size(J3)),J2];
%         elseif mode == 3
%             J = [Psi2.*J1./(Psi1.^2),J2./Psi1];
%         end
%     else
        if mode ==1 
            J=[-J1,J2];
        elseif mode == 2
            J=[J1,zeros(size(J1));zeros(size(J1)),J3];
        elseif mode == 3
            J = [J1,J2];
        end
%     end


    % Calculation of the desired prior-matrix with the functiom calc_RtR_prior 
    inv_model.jacobian_bkgnd.value=real(bkgnd1);
    RtR = calc_RtR_prior( inv_model );

    %The prior matrix is adapted to the formulation of the problem with
    %gamma_0^(tot) (cf Masterarbeit Seite 27)
    RtR=blkdiag(RtR,RtR);
    J_sym=J'*J;
    %JtJ results in mixed gradients of Jhigh and Jlow -> set to zero
    J_sym(1:nb_elem,(nb_elem+1):2*nb_elem) = 0;
    J_sym((nb_elem+1):2*nb_elem,1:nb_elem) = 0;

    % Definition of the reconstruction matrix
    RM=(J_sym + hp*RtR)\J';
end