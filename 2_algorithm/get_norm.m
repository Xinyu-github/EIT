function [residual,regularizer] = get_norm(imgr)
%% Residual and Regularization Norm is calculated here
    for i = 1:length(imgr)
        img = imgr(i);
        nb_elem = size(img.elem_data,1);
        residual_diff(i) = (norm(img.Psi2w-img.Psi1 -img.data2+img.alpha.*img.data1))/2;
        residual_abs(i) = img.beta*(norm(img.Psi1-img.data1)+norm(img.Psi2-img.data2))/2; %
        residual(i) = residual_diff(i)+residual_abs(i);
        L_diff = img.L(1:nb_elem,1:nb_elem)*img.elem_data/2;
        regularizer(i) = norm(L_diff);
        labels(i) = img.hp;
    end
        
    % Res = residual(1):(residual(end)-residual(1))/100:residual(end);
    % Reg = interp1(residual,regularizer,Res);

    % xx=log10(Res);
    % df1=gradient(Reg,xx);
    % df2=gradient(df1,xx);
    % id=sign(df2);
    % idx=strfind(id,[1 -1]);

    
%         figure
%         loglog(residual,regularizer);
%         xlabel('complete residual');
%         ylabel('regularizer');
%         for j = 1:length(imgr)
%         % text(residual,regularizer,num2str(labels));
%             txt = strcat('\lambda = ',num2str(labels(j)));   
%             text(residual(j),regularizer(j),txt,'HorizontalAlignment','right');
%         end
%         figure
%         loglog(residual_abs,regularizer);
%         % loglog(Res,Reg);
%         xlabel('residual abs');
%         ylabel('regularizer');
%         for j = 1:length(imgr)
%         % text(residual,regularizer,num2str(labels));
%             txt = strcat('\lambda = ',num2str(labels(j)));   
%             text(residual(j),regularizer(j),txt,'HorizontalAlignment','right');
%         end
%         figure
%         loglog(residual_diff,regularizer);
%         % loglog(Res,Reg);
%         xlabel('residual diff');
%         ylabel('regularizer');
%         for j = 1:length(imgr)
%         % text(residual,regularizer,num2str(labels));
%             txt = strcat('\lambda = ',num2str(labels(j)));   
%             text(residual(j),regularizer(j),txt,'HorizontalAlignment','right');
%         end
%         figure
%         loglog(residual_diff,regularizer+residual_abs);
%         % loglog(Res,Reg);
%         xlabel('residual diff');
%         ylabel('regularizer + residual abs');
%         for j = 1:length(imgr)
%         % text(residual,regularizer,num2str(labels));
%             txt = strcat('\lambda = ',num2str(labels(j)));   
%             text(residual(j),regularizer(j),txt,'HorizontalAlignment','right');
%         end

end
