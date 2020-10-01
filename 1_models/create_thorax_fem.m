function fmdl = create_thorax_fem(n_elec, nSkip, max_grid_size, grid, height)
% This function creates a thorax model with the given parameters
% typical values:
% n_elec = 32;
% nksip = 3;
% max_grid_size = 0.03;
% grid specifies which organ shapes should be included in the extruded model:
%   0 - only thorax contour, no organ shapes
%   1 - rlung, llung
%   2 - rlung
%   3 - llung
% height = 0;

%Loading Thorax , left and right lung from the library
thorax = shape_library('get','adult_male','boundary');
rlung  = shape_library('get','adult_male','right_lung');
llung  = shape_library('get','adult_male','left_lung');

if grid == 0
    shape = { height,{thorax},[4,19],max_grid_size};
elseif grid == 1
    shape = { height,{thorax,rlung,llung},[4,19],max_grid_size};
elseif grid == 2
    shape = { height,{thorax,rlung},[4,19],max_grid_size};    
elseif grid == 3
    shape = { height,{thorax,llung},[4,19],max_grid_size};
% elseif grid == 5
%     shape = { height,{thorax},[4,19],max_grid_size};
%     extra = { 'ball', strcat('solid ball = sphere(',num2str(cyl_center(1)),',0,0.5;0.1);')}; %radius 0.05*2*1 for thorax x-radius
%     %extra={'twoball','solid twoball = sphere(0.2,0.2,.7;0.2) or sphere(-0.2,-0.2,.4;0.3);'};
end

% Electrode Positions and Shape for 2D and 3D Stimulation
elec_pos = [ n_elec,1,height/2];
elec_shape = [0.05,0, 1 ];

%Compile of the forward Model with Netgen
if grid == 5
    [fmdl, mat_idx] = ng_mk_extruded_model(shape, elec_pos, elec_shape, extra);
else
    [fmdl, mat_idx] = ng_mk_extruded_model(shape, elec_pos, elec_shape);
end

[stim,msel] = mk_stim_patterns(n_elec,1,[0,nSkip+1],[0,nSkip+1],{'no_meas_current'},1);
fmdl.stimulation = stim; 
fmdl.mat_idx = mat_idx;

end