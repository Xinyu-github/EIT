function func_selected = my_create_simple_organ_shape(heart_center,heart_radius, spine_center)

% heart_radius = 0.15;%0.35; % Radius of the heart (modeled by a circle)
lung_limit_x = 0.10; % Value of x where the lung stops (the lung is modelled by a cropped circle), the lungs only spans for values of x>lung_limit_x
lung_limit_y_left = -0.7; % Parameter useful to model the hemothorax: the left lung spans only for values of y which are superior to lung_limit_y_left. The rest of what should be the lung is filled with blood.
lung_limit_y_right = -0.7; % Parameter useful to model the hemothorax: the right lung spans only for values of y which are superior to lung_limit_y_left. The rest of what should be the lung is filled with blood.
lung_radius = 0.35; % Radius of the circles that make the lungs
spine_radius_y = 0.1; % The spine is an ellipse, it therefore has a radius in y 
spine_radius_x = 0.1; % and a radius in x.
artefact_center = [-0.3,0]; artefact_radius = 0.2;
% spine_radius_z = 0.3;

% Selection function defining the domains where the organs span.
func_selected.heart = @(x,y,z) (x-heart_center(1)).^2+(y-heart_center(2)).^2 < heart_radius^2;

func_selected.spine = @(x,y,z) ((x-spine_center(1))/spine_radius_x).^2+((y-spine_center(2))/spine_radius_y).^2 <= 1.01 ;

% func_selected.lung = @(x,y,z) (  (((x-lung_center_right(1)).^2+(y-lung_center_right(2)).^2<lung_radius^2)&(x>lung_limit_x)&...
%                                   (y>lung_limit_y_right))  |  (((x+lung_center_right(1)).^2+(y-lung_center_right(2)).^2<lung_radius^2)&...
%                                   (x<-lung_limit_x)&(y>lung_limit_y_left)))&(~(func_selected.heart(x,y,z)))&(~(func_selected.spine(x,y,z)));
% func_selected.lung = @(x,y,z) (  (((x-lung_center_right(1)).^2+(5.*(y-lung_center_right(2))./4 - sqrt(abs(x-lung_center_right(1)))).^2<1)&(x>lung_limit_x)&...
%                                   (y>lung_limit_y_right))  |  (((x+lung_center_right(1)).^2+(5.*(y-lung_center_right(2))./4 - sqrt(abs(x+lung_center_right(1)))).^2<1)&...
%                                   (x<-lung_limit_x)&(y>lung_limit_y_left)))&(~(func_selected.heart(x,y,z)))&(~(func_selected.spine(x,y,z)));

% func_selected.blood = @(x,y,z) (((x+lung_center_right(1)).^2+(y-lung_center_right(2)).^2<lung_radius^2)&(x<-lung_limit_x)&...
%                                  (y<lung_limit_y_left)) & (~(func_selected.heart(x,y,z)))|(((x-lung_center_right(1)).^2+...
%                                  (y-lung_center_right(2)).^2<lung_radius^2)&(x>lung_limit_x)&(y<lung_limit_y_right))&...
%                                  (~(func_selected.heart(x,y,z)));
% func_selected.artefact = @(x,y,z) (x-artefact_center(1)).^2+(y-artefact_center(2)).^2 < artefact_radius^2;