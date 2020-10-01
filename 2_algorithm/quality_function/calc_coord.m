function coord=calc_coord(vect,npoints)
% Function to convert the coordinates of a matrix considered as a vector
% into their line/column equivalent

%ss57
coord=[floor(vect/npoints)+1,vect-npoints*(floor(vect/npoints))];

end
