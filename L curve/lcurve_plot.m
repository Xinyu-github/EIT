% load data
load('..\.\Reconstruction\recon_atright_beta0-0.1.mat');

% load('..\.\Reconstruction\recon_healthy_beta0-0.1.mat');

%% 2D plot
for jj = 1:length(Img)
[res,reg] = generate_lcurve(Img{jj});
end

%% 3D plot
Res = zeros(101,11);
Reg = zeros(101,11);
for jj = 1:length(Img)
    imgr = Img{jj};
    for i = 1:length(imgr)
        img = imgr(i);
        nb_elem = size(img.elem_data,1);
        residual(i) = norm(img.Psi2w-img.Psi1 -img.data2+img.alpha.*img.data1)+img.beta*(norm((img.Psi1-img.data1))+norm((img.Psi2-img.data2)));
        L_diff = img.L(1:nb_elem,1:nb_elem)*img.elem_data;
        regularizer(i) = norm(L_diff);
        labels(i) = img.hp;
    end
    residual = normalize(residual,'range');
    regularizer = normalize(regularizer,'range');

%     interpolate to smooth curve

    Res(:,jj) = residual(1):(residual(end)-residual(1))/100:residual(end);
    Reg(:,jj) = interp1(residual,regularizer,Res(:,jj));
end
B = repmat(Beta,101,1);
surf(Res,B,Reg)
xlabel('residual norm');
ylabel('Beta');
zlabel('regularisation norm');
title('3D Lcurve for healthy lung for hp \in [0.01,0.5]');
