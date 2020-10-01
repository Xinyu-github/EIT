function plot_specific_fdEITimage(pathology, beta_num, hp_num, save_fig)
    beta=[0.25,0.5,0.75,1,1.25,2.5,5,7.5,10];
    load(strcat('recon_c',num2str(pathology),'_B',num2str(beta_num),'.mat'));
    img = img_solve;
    img(hp_num).show_phase_contour_slice.abs_clim = 0.1;
    img(hp_num).show_phase_contour_slice.abs_ref = 0.1;
    %img(hp_num).show_phase_contour_slice.abs_clim = 100;
    %img(hp_num).show_phase_contour_slice.abs_ref = 50;
    img(hp_num).show_phase_contour_slice.full_phase = 0;
    img(hp_num).show_phase_contour_slice.phase_clim = [0, 90];
    img(hp_num).show_phase_contour_slice.label = 1;
    figure;
    set(gcf, 'Position', [100, 200, 650, 380],'Color','white');
    %img(hp_num).show_phase_contour_slice.fname = strcat(strrep(strcat('recon_c1_B',num2str(beta(beta_num)),'hp',num2str(hp_param(ii))),'.','-'),'.pdf');
    img(hp_num).show_phase_contour_slice.switch_contour =0;
    show_phase_contour_slice(img(hp_num),64,10);%title('1) Reference img');
    if save_fig
        savefig(strcat(strrep(strcat('recon_c',num2str(pathology),'_B',num2str(beta(beta_num)),'hp',num2str(hp_param(hp_num))),'.','-'),'.fig'));
        %Show single magnitude and phase image
        load 'colormaps.mat'
        fSize = 16;
        figure;
        img_abs=img(hp_num);
        img_abs.elem_data = abs(img_abs.elem_data);
        img_abs.calc_colours.cmap_type = cm_0.abs;
        show_slices(img_abs);
        cb1=eidors_colourbar(img_abs);
        ylabel(cb1,'Magnitude [S/m]');
        set(gca,'FontSize',fSize);
        savefig(strcat(strrep(strcat('recon_c',num2str(pathology),'_B',num2str(beta(beta_num)),'hp',num2str(hp_param(hp_num)),'_magnitude'),'.','-'),'.fig'));

        figure;
        img_phase=img(hp_num);
        img_phase.elem_data = angle(img_phase.elem_data)/pi*180;
        show_slices(img_phase);
        %include gray bg
        cm_phase_mod=cm_3.phase;
        cm_phase_mod(1,:)=[0.349 0.498 0.498];
        colormap(cm_phase_mod);
        img_phase.calc_colours.ref_level =  50;%  centre of the colour scale
        img_phase.calc_colours.clim      =  50;%  max diff from ref_level
        %img_abs.calc_colours.cmap_type = cm_0.phase;
        cb2=eidors_colourbar(img_phase);
        ylabel(cb2,'Phase [°]');
        set(gca,'FontSize',fSize);
        savefig(strcat(strrep(strcat('recon_c',num2str(pathology),'_B',num2str(beta(beta_num)),'hp',num2str(hp_param(hp_num)),'_phase'),'.','-'),'.fig'));
    end
end
