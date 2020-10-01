function center=calc_center(vect,npoints)
% Function to calculate the coordinates of the isobarycentre of the vector
% with the origin of the frame being the center of the "circular" grid

% ss58
coord_vect=calc_coord(vect,npoints); %ss57
center=[mean(coord_vect(:,1)-0.5),mean(coord_vect(:,2)-0.5)]-[npoints/2,npoints/2];

end
