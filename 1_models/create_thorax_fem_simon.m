function fmdl = create_thorax_fem_simon(n_elec, nSkip, max_grid_size, lung_grid, height, cyl_center)
% This is the same funcation as 'my_create_thorax_fem' stored in './models'
% and has been copied to this location and renamed in order to prevent
% filname and function name problems.

%% create thorax model
% n_elec = 32; max_grid_size = 0.03;
%Loading Thorax , left and right lung from the library
thorax = shape_library('get','adult_male','boundary');
rlung  = shape_library('get','adult_male','right_lung');
llung  = shape_library('get','adult_male','left_lung');
llung1  = shape_library('get','neonate','left_lung');

if lung_grid == 1
    shape = { height,{thorax,rlung,llung},[4,19],max_grid_size};
elseif lung_grid == 2
    shape = { height,{thorax,rlung},[4,19],max_grid_size};    
elseif lung_grid == 3
    shape = { height,{thorax,llung},[4,19],max_grid_size};
elseif lung_grid == 4
    shape = { height,{thorax,rlung,llung,[llung1(:,1)*0.6+0.25,llung1(:,2)*0.6]},[4,19],max_grid_size};
elseif lung_grid == 5
    shape = { height,{thorax},[4,19],max_grid_size};
    extra = { 'ball', strcat('solid ball = sphere(',num2str(cyl_center(1)),',0,0.5;0.1);')}; %radius 0.05*2*1 for thorax x-radius
    %extra={'twoball','solid twoball = sphere(0.2,0.2,.7;0.2) or sphere(-0.2,-0.2,.4;0.3);'};
else
    shape = { height,{thorax},[4,19],max_grid_size};
end

% Electrode Positions and Shape for 2D and 3D Stimulation
elec_pos = [ n_elec,1,height/2];
elec_shape = [0.05,0, 1 ];

%Compile of the forward Model with Netgen
if lung_grid == 5
    [fmdl, mat_idx] = ng_mk_extruded_model(shape, elec_pos, elec_shape, extra);
else
    [fmdl, mat_idx] = ng_mk_extruded_model(shape, elec_pos, elec_shape);
end

[stim,msel] = mk_stim_patterns(n_elec,1,[0,nSkip+1],[0,nSkip+1],{'no_meas_current'},1);
fmdl.stimulation = stim; 
fmdl.mat_idx = mat_idx;

end