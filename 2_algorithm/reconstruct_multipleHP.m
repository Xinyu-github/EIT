function [img_solve] = reconstruct_multipleHP(fmdl_recon, hp_param,beta,weighted,init_iter, v1, v2, mode,title)
% RECONSTRUCT_MULTIPLEHP: Reconstructs a set of hp_param vor multiple set
% of beta variables
    for k=1:length(beta(:,1))
        clear img_solve
        for ii = 1:length(hp_param)
            imdl    = set_weighted_invprob_properties(fmdl_recon, init_iter, hp_param(ii),...
                weighted,beta(k,ii),mode);
            img_solve(ii)= inv_solve(imdl, v1, v2);

            disp('o-----------------------------------------------------------------------------------o')
            disp(strcat('Finished_', title ,' reconstruction for hyperparameter: '));
            disp(strcat(num2str(ii),'/',num2str(length(hp_param)),' -> hp=',num2str(hp_param(ii))));
            disp('with beta: ')
            disp(strcat(num2str(k),'/',num2str(length(beta(:,1))),' -> beta=',num2str(beta(k,ii))));
            disp('o-----------------------------------------------------------------------------------o')
        end
        beta_param = beta(k,:);
        save(strcat('./3_reconstructions/results/recon_',title,'_allHP.mat'),'img_solve','beta_param','hp_param');
    end
end