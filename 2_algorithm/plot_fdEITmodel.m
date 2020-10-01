function plot_fdEITmodel(img_diff, title_, save_fig,name,fname)
%     img(hp_num).show_phase_contour_slice.abs_clim = 0.1;
%     img(hp_num).show_phase_contour_slice.abs_ref = 0.1;
%     %img(hp_num).show_phase_contour_slice.abs_clim = 100;
%     %img(hp_num).show_phase_contour_slice.abs_ref = 50;
%     img(hp_num).show_phase_contour_slice.full_phase = 0;
%     img(hp_num).show_phase_contour_slice.phase_clim = [0, 90];
%     img(hp_num).show_phase_contour_slice.label = 1;


%  ref=0.1;
%  clim = 0.1;
%  img_diff.show_phase_contour_slice.abs_clim = clim;
%  img_diff.show_phase_contour_slice.abs_ref = ref;
img_diff.show_phase_contour_slice.full_phase = 0;
img_diff.show_phase_contour_slice.phase_clim = [0, 90];
img_diff.show_phase_contour_slice.label = 1;

%Export figure
%img_diff.show_phase_contour_slice.fname = "reference_healthy_ellipsoid";
figure,set(gcf, 'Position', [100, 200, 650, 380],'Color','white');
show_phase_contour_slice(img_diff,64,10);%title('1) Reference img');
title(title_);
    if save_fig
        
        %savefig(strcat(name,'.fig'));
        saveas(gcf,fullfile(fname, name),'jpg');
        return;
        %Show single magnitude and phase image
        load 'colormaps.mat'
        fSize = 16;
        figure;
        img_abs=img_diff;
        img_abs.elem_data = abs(img_abs.elem_data);
        img_abs.calc_colours.cmap_type = cm_0.abs;
        img_abs.calc_colours.clim = clim;
        img_abs.calc_colours.ref_level = ref;
        show_slices(img_abs);
        title([title_,'mag']);
        cb1=eidors_colourbar(img_abs);
        ylabel(cb1,'Magnitude [S/m]');
        set(gca,'FontSize',fSize);
        savefig(strcat(name,'_magnitude.fig'));
        
        figure;
        img_phase=img_diff;
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
        title([title_,'phase']);
        savefig(strcat(name,'_phase.fig'));
    end
end
