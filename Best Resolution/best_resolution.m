nSkip = 0;
load('all_cases.mat');
fmdl_recon = create_thorax_fem_simon(32, nSkip, 0.05, 0, 0, 0);
hp = logspace(-7,-5,10);
beta = 0.5;
BR = zeros(length(beta),length(hp));
for i = 1:length(beta)
    for j = 1:length(hp)
        imdl    = mk_hpselection_imdl(cases{i}.fmdl,fmdl_recon, 1,0,'bestres',1,0,2);
        BR(i,j) = calc_image_br(imdl,hp(j));
    end
end