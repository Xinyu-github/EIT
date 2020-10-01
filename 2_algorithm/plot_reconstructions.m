function plot_reconstructions(img_solve,iter,name)
%PLOT_RECONSTRUCTIONS: plots all given images specified by iter
%   Detailed explanation goes here
    for i = iter
        figure;set(gcf, 'Position', [100, 200, 650, 380],'Color','white');
        string_hp = num2str(img_solve(i).hp);
        string_beta = num2str(img_solve(i).beta);
        show_phase_contour_slice(img_solve(i),64,10);
        title(strcat(name,' - hp: ',string_hp,' - beta: ',string_beta));
    end
end

