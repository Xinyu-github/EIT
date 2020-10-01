function displacement=calc_displacement(detected_organ,center_ref,npoints, radius_pix)
% Function to compute the displacement for a detected organ

%ss61
    center_organ=calc_center(detected_organ,npoints);
    displacement=norm(center_organ-center_ref)/radius_pix;
    %Should the displacement be truncated to 1 if it is greater than 1 ?
end
