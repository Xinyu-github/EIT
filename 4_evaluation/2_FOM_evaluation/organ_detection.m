function [BWmatrix,detection_infos,selected_BW,selected_indx]=organ_detection(recons_matrix,min_wind,max_wind,organ_grid)
% Function to detect the organs i.e. detect the pixels of a matrix whose
% value is inside a given range and detect the shape that the pixels make on
% the domain (i.e. the image of the organs)

%ss56
% recons_cond within [min_cond_wind, max_cond_wind]
    BWmatrix=recons_matrix>=min_wind & recons_matrix<=max_wind; %recons_matrix pixels locate within moving window
    detection_infos=bwconncomp(BWmatrix,8);  %find connecting nbh(target) on the selected recons_matrix pixels  
% recons_region intersects ref_organ
    [selected_BW,selected_indx.recons,selected_indx.ref] = intersect(recons_matrix,organ_grid,'stable');
end
