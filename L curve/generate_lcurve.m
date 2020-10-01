function [residual,regularizer] = generate_lcurve(imgr)

for i = 1:length(imgr)
    img = imgr(i);
    nb_elem = size(img.elem_data,1);
    if  ~isfield(img,'mode')||img.mode ==1 
        residual(i) = norm(img.Psi2w-img.data2-img.Psi1 +img.alpha*img.data1)+img.beta*(norm(img.Psi2-img.data2)+norm(-img.Psi1 +img.data1));
        
%         residual(i) = norm((img.Psi2w-img.Psi1 )./(img.data2-img.alpha*img.data1));
    elseif img.mode == 2
       residual(i) = norm(img.Psi2-img.data2)+norm(img.alpha*(img.Psi1 -img.data1));
%         residual(i) = norm((img.Psi2w+img.Psi1 -img.data2-img.alpha*img.data1)./(img.data2+img.alpha*img.data1));
%         residual(i) =  norm(img.Psi2w-img.Psi1 -img.data2+img.alpha*img.data1);
    elseif img.mode==3
        residual(i) = norm(img.Psi2w-img.data2-img.Psi1 +img.alpha*img.data1);
    end
    cond = [img.cond1;img.cond2];
     regularizer(i) = norm(img.L*cond);
%      regularizer(i) =  norm(img.Psi2-img.data2)+norm(img.Psi1 -img.data1);
   
    labels(i) = img.hp;
end
% Res = residual(1):(residual(end)-residual(1))/100:residual(end);
% Reg = interp1(residual,regularizer,Res);

% xx=log10(Res);
% df1=gradient(Reg,xx);
% df2=gradient(df1,xx);
% id=sign(df2);
% idx=strfind(id,[1 -1]);

figure
loglog(residual,regularizer);
% loglog(Res,Reg);
xlabel('residual');
ylabel('regularizer');
% for j = 1:length(imgr)
% % text(residual,regularizer,num2str(labels));
%     txt = strcat('\lambda = ',num2str(labels(j)));   
%     text(residual(j),regularizer(j),txt,'HorizontalAlignment','right');
% end

for j = 1:length(imgr)
% text(residual,regularizer,num2str(labels));

    txt = strcat('\lambda = ',num2str(labels(j)));   
    text(residual(j),regularizer(j),txt,'HorizontalAlignment','right');

end
