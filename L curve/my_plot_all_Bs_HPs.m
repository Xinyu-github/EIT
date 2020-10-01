Beta=0:0.01:0.1;%1];%,1.25,2.5,5,7.5,10];
init_val=0.1+0.01i;
% load('recon_healthy_beta0-0.1.mat');
% load('recon_atright_beta0-0.1.mat');
for jj = 2:3
%     load(strcat('recon_c1_B',num2str(jj),'.mat'));
hp_param = Hp;
img = img_solve(12:27);
%norm_noise = norm(periodogram(([1 1j]*randn(2,length([img(i).cond1; img(i).cond2]))/sqrt(2))'));
for i=1:length(hp_param)%number of hyper-parameter
%     img_fwd_model = mk_image(img(i).inv_model);
%     img_fwd_model.elem_data=img(i).alpha.*img(i).cond2;
%     Psi2w=fwd_solve(img_fwd_model);
    Psi2w = img(i).Psi2w;
    ls_residual(i) = norm((Psi2w-(img(i).Psi1))-(img(i).data2-(img(i).alpha.*img(i).data1)));%/norm((img(i).alpha.*img(i).data1));
    ls_mixed(i) = (norm((img(i).Psi1-img(i).data1))+norm((img(i).Psi2-img(i).data2)));%/norm(img(i).data1));  
    ls_term(i) = ls_residual(i)+img(i).beta*ls_mixed(i);
    %[pxx,w] = periodogram((Psi2w.meas-(img(i).Psi1))-(img(i).data2-(img(i).alpha.*img(i).data1)));
    %[pxx2,w] = periodogram(img(i).Psi1-img(i).data1);
    %[pxx3,w] = periodogram(img(i).Psi2-img(i).data2);
    %norm_sol(i) = norm(pxx)+norm(pxx2)+norm(pxx3);
    %figure;
    %plot(w,10*log10(pxx));
    reg_term(i) = norm(img(i).L*([img(i).cond1; img(i).cond2]));%-init_val));%/norm(cond{1,i}(1:length(cond{1,i})/2));
    %reg_term(i) = norm(img(i).cond1)+norm(img(i).cond2);
    %figure(400+i),set(gcf, 'Position', [100, 200, 650, 380],'Color','white');
    %show_phase_contour_slice(img(i),64,10);%title('1) Reference img');
end
%figure(200)
%plot(hp_param, norm_sol,'-s','linewidth',2,'DisplayName',strcat('\beta=',num2str(beta(jj)))),hold on,

%figure(201)
%loglog(ls_mixed(:), reg_term(:),'-s','linewidth',2,'DisplayName',strcat('\beta=',num2str(beta(jj)))),hold on,

%% Calculation of Curvature for all selected Bs
HPindx_maxK = []; % collects max curvature of each L-curve (at corner)
dx  = gradient(ls_term(1:end));
ddx = gradient(dx);
dy  = gradient(reg_term(1:end));
ddy = gradient(dy);
num   = abs(dx .* ddy - ddx .* dy);
denom = dx .* dx + dy .* dy;
denom = sqrt(denom);
denom = denom.^3; %* denom * denom;
curvature(jj,:) = num ./ denom;
[val,indx] = max(curvature(jj,(1:end-1)));
opt_hp_param(jj) = hp_param(indx-1);

figure(202);
loglog(ls_term(:), reg_term(:),'-s','linewidth',2,'DisplayName',strcat('\beta=',num2str(Beta(jj)))),hold on,
figure(203);
plot(hp_param(1:end), curvature(jj,:), 'Linewidth', 2,'DisplayName',strcat('\beta=',num2str(Beta(jj)))), hold on,
display(strcat('beta=',num2str(Beta(jj)),' max curve at lambda=',num2str(opt_hp_param(jj))));
% figure(jj),set(gcf, 'Position', [100, 200, 650, 380],'Color','white');
% show_phase_contour_slice(img(indx-1),64,10);%title('1) Reference img');
end

figure(204);
legend('show')
figure(205);
legend('show')


% beta=0.25 max curve at lambda=
% beta=0.5 max curve at lambda=0.9 i=8
% beta=0.75 max curve at lambda=
%choose beta=0.5 max curve at lambda=0.9 i=8
% plot_specific_fdEITimage(1, 1, 4, 1) %corner point1
