function img = beta_selection(img_diff,img_abs,PlotAll)

% Input:
%     img_diff: the best reconstructed image using differential norm
%     img_abs: the best reconstructed image using abs norm
%
% Output:
%     img:  image that the conductivity is generated with best beta


    beta = 0:0.1:1;
    img = img_diff;
    for i = 1:length(beta)
        cond1 = (img_diff.cond1 + beta(i)*img_abs.cond1)/(1+beta(i));
        img.elem_data = cond1;
        Psi1 = fwd_solve(img);
        Psi1 = Psi1.meas;
        cond2 = (img_diff.cond2 + beta(i)*img_abs.cond2)/(1+beta(i));
        img.elem_data = img.alpha*cond2;
        Psi2w = fwd_solve(img);
        Psi2w = Psi2w.meas;
        img.elem_data = cond2;
        Psi2 = fwd_solve(img);
        Psi2 = Psi2.meas;
        residual(i) = norm((Psi2w-Psi1-img.data2+img.alpha*img.data1));
        residual2(i) = (norm(Psi2-img.data2)+norm(img.alpha*(Psi1-img.data1)));
    %     residual2(i) = residual2(i)/(1+beta(i));
    end
    img.res1 = residual;
    img.res2 = residual2;
    r = sqrt((residual/max(residual)).^2+(residual2/max(residual2)).^2);
    i = find(r==min(r));
    img.elem_data = (img_diff.cond2 + beta(i)*img_abs.cond2)/(1+beta(i))-(img_diff.cond1 + beta(i)*img_abs.cond1)/(1+beta(i));
    figure;
        plot(beta,r);
%         i = find(residual==min(residual));
    
    if PlotAll
        for i = 1:length(beta)
            cond1 = (img_diff.cond1 + beta(i)*img_abs.cond1)/(1+beta(i));
            cond2 = (img_diff.cond2 + beta(i)*img_abs.cond2)/(1+beta(i));
            img.elem_data = (img_diff.cond2 + beta(i)*img_abs.cond2)/(1+beta(i))-(img_diff.cond1 + beta(i)*img_abs.cond1)/(1+beta(i));
            plot_fdEITmodel(img,strcat('beta=',num2str(beta(i))),0);
        end
%         
    
        
        
    end


end