function fmdl = my_create_thorax_fem(n_elec,max_grid_size)
%% create thorax model
% n_elec = 32; max_grid_size = 0.03;
%Loading Thorax , left and right lung from the library
thorax = shape_library('get','adult_male','boundary');
rlung  = shape_library('get','adult_male','right_lung');
llung  = shape_library('get','adult_male','left_lung');

shape = { 0,{thorax+0.5,rlung+0.5,llung+0.5},[4,19],max_grid_size};

% Electrode Positions and Shape for 2D and 3D Stimulation
elec_pos = [ n_elec, 0];
elec_shape = [0.05, 0, 1]; 

%Compile of the forward Model with Netgen
[fmdl, mat_idx] = ng_mk_extruded_model(shape, elec_pos, elec_shape);
[stim,msel] = mk_stim_patterns(n_elec,1,[0,4],[0,4],{'no_meas_current'},1);
fmdl.stimulation = stim; 
fmdl.mat_idx = mat_idx;

end