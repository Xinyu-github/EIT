function show_phase_contour_slice( img, res, kont_num , varargin)
%% line47
%show_phase_contour_slice( img,res,kont_num, z_cuts, x_cuts, y_cuts)
%  works like show_3d_slices just for contourslice()
%  z_cuts = planes in z to do a cut
%  x_cuts = planes in x to do a cut
%  y_cuts = planes in y to do a cut
% Default show 2 z_cuts and 1 x and 1 y cut
%img - Eidors img-object
%res - resolution of the contour plot
%kont_num - number of contour-lines
%kont_num can be a vector of phasevalues where contour-slices should be
%drawn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Options:
%img.show_phase_contour_slice.full_phase:
%if true then the whole phase is displayed from -180 up to 180
%
%img.show_phase_contour_slice.abs_cm and
%img.show_phase_contour_slice.phase_cm:
%use this to determine the colormaps of the absolute and phase value. predetermined
%colormaps are in the file "colormaps.mat". the default colormap is cm_1
%but everycolormap in matlab can be used
%
%img.show_phase_contour_slice.line_width: 
%determine the thickness of the contourplot of the phasevalue (default = 2)
%
%These Options only work with 2D-Models. 3D will be added soon!


% img, res and kont_num are the only things to set when working with
% 2D-models

fSize = 16;

try
    img.show_phase_contour_slice.switch_contour;
catch
    img.show_phase_contour_slice.switch_contour = 0;
end

elem_data = img.elem_data;
if img.show_phase_contour_slice.switch_contour    
    elem_data_phase = abs(img.elem_data);               % absolute value
    elem_data_abs   =radtodeg(angle((img.elem_data)));  % phase in degrees
else
    elem_data_phase = radtodeg(angle((img.elem_data))); % phase in degrees
    elem_data_abs   =abs(img.elem_data);                % absolute value 
end

% try %Use max_abs_ref to scale the absolute values to the reference maximum
%     img.show_phase_contour_slice.max_abs_ref; %max_abs_ref=max(abs(img_diff.elem_data))
%     elem_data_abs = elem_data_abs-min(elem_data_abs);
%     elem_data_abs = elem_data_abs/max(elem_data_abs);
% 
%     elem_data_abs = elem_data_abs.*img.show_phase_contour_slice.max_abs_ref;
%     %elem_data_abs = elem_data_abs+min(abs_ref);
% catch
% end
img.elem_data = elem_data_abs;

try
    img.show_phase_contour_slice.label;
catch
    img.show_phase_contour_slice.label = 1;
end

try
    img.show_phase_contour_slice.full_phase;
catch
    img.show_phase_contour_slice.full_phase = 0;
end

load 'colormaps.mat'
try
    img.show_phase_contour_slice.abs_cm;
    img.show_phase_contour_slice.phase_cm;
catch
    img.show_phase_contour_slice.abs_cm = cm_0.abs;
    img.show_phase_contour_slice.phase_cm= cm_3.phase;
end

try 
    img.show_phase_contour_slice.line_width;
catch
    img.show_phase_contour_slice.line_width=2;
end

img.calc_colours.cmap_type = img.show_phase_contour_slice.abs_cm;

%%%%%%%%%%%%%%Credit to Sebastien Dambrun%%%%%%%%%%%%%%%%
%maybe comment it to set different ref-levels etc.
    [jnk,ref_lev,max_scale] = scale_for_display( img.elem_data);
%    max_scale = abs(max(elem_data_abs)-ref_lev);
%    ref_lev = 0.1271; %median(elem_data_abs)
%    ref_lev = median(elem_data_abs);
try 
    img.calc_colours.ref_level; 
catch
    img.calc_colours.ref_level = ref_lev;
end
try
    img.calc_colours.clim;
catch
    img.calc_colours.clim = max_scale;
end
try 
    np = img.calc_colours.npoints;
catch
    np = calc_colours('npoints');
end

mdl_min = min(img.fwd_model.nodes);%get min and max coordinates
mdl_max = max(img.fwd_model.nodes);
mdl_rng = mdl_max - mdl_min;
np = round(mdl_rng/min(mdl_rng) * np);
%%%%%%%%%%%%%%Credit to Sebastien Dambrun%%%%%%%%%%%%%%%%

fwd_model_dim = size(img.fwd_model.nodes);

if fwd_model_dim(1,2) ==3
    [x_cuts, y_cuts,z_cuts] = get_cuts(img,varargin{:});%get cuts for the mdl_slices mapper with a subfunction



    img.elem_data = elem_data_abs;
    %%


    ax2 = gca;              %use second coordinate system
    %%%%%%%%%%%%%%Credit to Sebastien Dambrun%%%%%%%%%%%%%%%%


    show_3d_slices(img, z_cuts,x_cuts, y_cuts);   %Draw the abs as a normal slice
    freezeColors                                  %Freeze the colors of the already drawn plots credit to John Iversen

    ax1 = axes;
    cb1 = eidors_colourbar(img);                  %create the colorbar for the absolute values
    set(cb1,'Position',[.1 .11 .0250 .815]);      %colorbar to the left

    ylabel(cb1,'Absolute Value');                 
    axis off                                      %the second axes are just to get the second Colorbar
                                                  %everything else is turned
                                                  %off

    axes(ax2);
    if fullphase == 0
    set(ax2, 'CLim', [min(elem_data_phase), max(elem_data_phase)]); %Set min and max limitations
    else
    set(gca, 'CLim', [-180, 180]); %Set min and max limitations if desired
    end
    colormap(ax2,'cool');        %use second color palette
    cb2 = colorbar;              %create second colorbar on the right for phase values
    ylabel(cb2,'Phase Value');

    img.elem_data = elem_data_phase;
    img.fwd_model.boundary = find_boundary(img.fwd_model);

    for i= 1:length(z_cuts)
        %Prepare mdl_slice_mapper function to return all necessary nodes

    img.fwd_model.mdl_slice_mapper.level = [inf inf z_cuts(i)];
    img.fwd_model.mdl_slice_mapper.x_pts = linspace(mdl_min(1),mdl_max(1),res);
    img.fwd_model.mdl_slice_mapper.y_pts = linspace(mdl_min(2),mdl_max(2),res);

    q= mdl_slice_mapper(img.fwd_model,'elem'); %return nodes of a specific z-slice

    %Build the coordinate grid for the contourslice-function
    [x y z] =meshgrid(linspace(mdl_min(1),mdl_max(1),res),...
                        linspace(mdl_min(2),mdl_max(2),res),...
                            linspace(mdl_min(3),mdl_max(3),res));

    %%Have to use this so you can read in zeros values of elem_data
    for m=1:1:res
        for n=1:1:res
            for p=1:1:res
                if q(m,n)>0
                    z_data(m,n,p) = img.elem_data(q(m,n));
                else
                    z_data(m,n,p) = 0;
                end
            end
        end
    end
    contourslice(x,y,z,z_data,[],[],z_cuts(i),kont_num); %use the contourslice-plot to draw the phase values onto the absolute value slice
        hold on
    end

    for i= 1:length(x_cuts)
        %Prepare mdl_slice_mapper function to return all necessary nodes
    img.fwd_model.mdl_slice_mapper.level = [x_cuts(i) inf inf];
    img.fwd_model.mdl_slice_mapper.x_pts = linspace(mdl_min(2),mdl_max(2),res);
    img.fwd_model.mdl_slice_mapper.y_pts = linspace(mdl_min(3),mdl_max(3),res);

    r= mdl_slice_mapper(img.fwd_model,'elem'); %return nodes of a specific slice

    x_size = size(img.fwd_model.mdl_slice_mapper.x_pts,2);
    y_size = size(img.fwd_model.mdl_slice_mapper.y_pts,2);

    %Build the coordinate grid for the contourslice-function


    [x y z] =meshgrid( linspace(mdl_min(1),mdl_max(1),x_size),...
                        linspace(mdl_min(2),mdl_max(2),y_size),...
                        linspace(mdl_min(3),mdl_max(3),res));

    %%Have to use this so you can read in zeros values of elem_data
    for m=1:1:res
        for n=1:1:res
            for p=1:1:res
                if r(p,n)>0
                    x_data(m,n,p) = img.elem_data(r(p,n));
                else
                    x_data(m,n,p) = 0;
                end
            end
        end
    end
    x_data = rot90(x_data,3); % rotate matrix to fit it right

    contourslice(x,y,z,x_data,x_cuts(i),[],[],kont_num); %use the contourslice-plot to draw the phase values onto the absolute value slice

        hold on
    end

    for i= 1:length(y_cuts)

    img.fwd_model.mdl_slice_mapper.level = [inf y_cuts(i) inf];
    img.fwd_model.mdl_slice_mapper.x_pts = linspace(mdl_min(1),mdl_max(1),res);
    img.fwd_model.mdl_slice_mapper.y_pts = linspace(mdl_min(3),mdl_max(3),res);

    s= mdl_slice_mapper(img.fwd_model,'elem'); %return nodes of a specific z-slice
    %Build the coordinate grid for the contourslice-function
    x_size = size(img.fwd_model.mdl_slice_mapper.x_pts,2);
    y_size = size(img.fwd_model.mdl_slice_mapper.y_pts,2);

    [x y z] =meshgrid(linspace(mdl_min(1),mdl_max(1),x_size),...
                        linspace(mdl_min(2),mdl_max(2),y_size),...
                            linspace(mdl_min(3),mdl_max(3),res));

    %%Have to use this so you can read in zeros values of elem_data
    for m=1:1:x_size
        for n=1:1:res
            for p=1:1:y_size
                if s(p,m)>0
                    y_data(m,n,p) = img.elem_data(s(p,m));
                else
                    y_data(m,n,p) = 0;
                end
            end
        end
    end
    y_data = rot90(y_data,1);% rotate matrix to fit it right
    contourslice(x,y,z,y_data,[],y_cuts(i),[],kont_num); %use the contourslice-plot to draw the phase values onto the absolute value slice

        hold on
    end

else

    img.fwd_model.mdl_slice_mapper.x_pts = linspace(mdl_min(1),mdl_max(1),res);
    img.fwd_model.mdl_slice_mapper.y_pts = linspace(mdl_max(2),mdl_min(2),res);

    x =linspace(mdl_min(1),mdl_max(1),res);
    y =linspace(mdl_min(2),mdl_max(2),res);


    q= mdl_slice_mapper(img.fwd_model,'elem'); %return nodes of a specific z-slice
    img.fwd_model.mdl_slice_mapper.x_pts = linspace(mdl_min(1),mdl_max(1),res);
    img.fwd_model.mdl_slice_mapper.y_pts = linspace(mdl_min(2),mdl_max(2),res);

    ax2 = gca;
    
    try
%    img.calc_colours.ref_level =  %  centre of the colour scale
    img.calc_colours.clim  = img.show_phase_contour_slice.abs_clim;
    img.calc_colours.ref_level  = img.show_phase_contour_slice.abs_ref;
%    img.eidors_colourbar.tick_vals = [-20:20]/10;
    end
    
    %show_fem(img);
    rimg = calc_slices(img); %Draw slice of the 2D model
    rimg = calc_colours(rimg,img);
    image([mdl_min(1),mdl_max(1)],[mdl_max(2),mdl_min(2)],rimg);
    
    %axis square
    pbaspect([mdl_max(1)-mdl_min(1) mdl_max(2)-mdl_min(2) 1])
    pb = (mdl_max(1)-mdl_min(1))/(mdl_max(2)-mdl_min(2));

    freezeColors
    
    pos = get(ax2,'Position');
    
    ax1 = axes;
    %set(ax1, 'CLim', [min(elem_data_abs), max(elem_data_abs]);
    
    cb1 = eidors_colourbar(img); %create the colorbar for the absolute values

    if img.show_phase_contour_slice.label == 1
        ylabel(cb1,'Magnitude [S/m]');        
    end
    axis off 

    axes(ax2);
    try
        set(ax2, 'CLim', img.show_phase_contour_slice.phase_clim)
    catch
        if img.show_phase_contour_slice.full_phase == 0
        set(ax2, 'CLim', [min(elem_data_phase), max(elem_data_phase)]); %Set min and max limitations
        %set(ax2, 'CLim', [0,100]); %Set min and max limitations
        else
        set(gca, 'CLim', [-180, 180]); %Set min and max limitations if desired
        end
    end
    
    colormap(ax2,img.show_phase_contour_slice.phase_cm);        %use second color palette
    %colormap(ax2,'cool');        %use second color palette
    cb2 = colorbar;              %create second colorbar on the right for phase values
    if round(pb,2) == 1
        set(cb1,'Position',[pos(1)+0.057 pos(2) 0.025 pos(4)]);%[.187 .11 .0250 .817]);   %colorbar to the left
        set(cb2,'Position',[pos(3)+0.048 pos(2) 0.025 pos(4)]);%set(cb2,'Position',[.823 .11 .0250 .817]);   %colorbar to the right
    elseif round(pb,2) == 1.47
        %set(cb1,'Position',[pos(1)-0.023 pos(2)+0.054 0.025 pos(4)-0.108]);%[.187 .11 .0250 .817]);   %colorbar to the left
        %set(cb2,'Position',[pos(4)+0.083 pos(2)+0.054 0.025 pos(4)-0.108]);%set(cb2,'Position',[.823 .11 .0250 .817]);   %colorbar to the right
        set(cb1,'Position',[pos(1)+0.014 pos(2)+0.0 0.025 pos(4)-0]);%[.187 .11 .0250 .817]);   %colorbar to the left
        set(cb2,'Position',[pos(4)+0.051 pos(2)+0.0 0.025 pos(4)-0]);%set(cb2,'Position',[.823 .11 .0250 .817]);   %colorbar to the right
    elseif round(pb,2) == 1.46
        %set(cb1,'Position',[pos(1)-0.023 pos(2)+0.054 0.025 pos(4)-0.108]);%[.187 .11 .0250 .817]);   %colorbar to the left
        %set(cb2,'Position',[pos(4)+0.083 pos(2)+0.054 0.025 pos(4)-0.108]);%set(cb2,'Position',[.823 .11 .0250 .817]);   %colorbar to the right
        set(cb1,'Position',[pos(1)+0.014 pos(2)+0.0 0.025 pos(4)-0]);%[.187 .11 .0250 .817]);   %colorbar to the left
        set(cb2,'Position',[pos(4)+0.051 pos(2)+0.0 0.025 pos(4)-0]);%set(cb2,'Position',[.823 .11 .0250 .817]);   %colorbar to the right
    else
        set(cb1,'Position',[pos(1)-0.02 pos(2) 0.025 pos(4)]);%[.187 .11 .0250 .817]);   %colorbar to the left
        set(cb2,'Position',[pos(4)+0.083 pos(2) 0.025 pos(4)]);%set(cb2,'Position',[.823 .11 .0250 .817]);   %colorbar to the right
    end
    if img.show_phase_contour_slice.label == 1
        ylabel(cb2,'Phase [°]');
    end
    
    

    img.elem_data = elem_data_phase;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not used afterwards, due to contourslice function
%      ref_lev = median(elem_data_phase);
% try 
%     img.calc_colours.ref_level; 
% catch
%     img.calc_colours.ref_level = ref_lev;
% end
% try
%     img.calc_colours.clim;
% catch
%     img.calc_colours.clim = ref_lev;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Build the coordinate grid for the contourslice-function
    [x y z] =meshgrid(linspace(mdl_min(1),mdl_max(1),res),...
                        linspace(mdl_min(2),mdl_max(2),res),...
                            linspace(mdl_min(2),mdl_max(2),res));

    %Have to use this so you can readin already pruned models
    for m=1:1:res
        for n=1:1:res
            for p=1:1:res
                if q(m,n)>0
                    z_data(m,n,p) = img.elem_data(q(m,n));
                else
                    z_data(m,n,p) = 0;
                end
            end
        end
    end
    c=contourslice(x,y,z,z_data,[],[],0,kont_num);
    set(c,'LineWidth',img.show_phase_contour_slice.line_width);
    %show_electrodes_2d(img.fwd_model, 1);
end
hold off
axis off
set(ax2,'FontSize',fSize);
set(ax1,'FontSize',fSize);

try
    export_fig(img.show_phase_contour_slice.fname, '-pdf');
catch
end

end

%Helpfunction from show_3d_slices
function [x_cuts, y_cuts, z_cuts] =  get_cuts(img, varargin)
   mdl_max= max(img.fwd_model.nodes);
   mdl_min= min(img.fwd_model.nodes);
   if nargin==1;
      % Default show 2 z_cuts and 1 x and 1 y cut
       x_cuts= linspace(mdl_min(1), mdl_max(1), 3); x_cuts([1,3])=[];
       y_cuts= linspace(mdl_min(2), mdl_max(2), 3); y_cuts([1,3])=[];
       z_cuts= linspace(mdl_min(3), mdl_max(3), 4); z_cuts([1,4])=[];
   elseif nargin==2;
       z_cuts= varargin{1};
       x_cuts= [];
       y_cuts= [];
   elseif nargin==3;
       z_cuts= varargin{1};
       x_cuts= varargin{2};
       y_cuts= [];
   elseif nargin==4;
       z_cuts= varargin{1};
       x_cuts= varargin{2};
       y_cuts= varargin{3};
   else 
       error('too many inputs');
   end 
end

function show_electrodes_2d(mdl, number_electrodes)
    if ~isfield(mdl,'electrode'); return; end

    ee= get_boundary( mdl );
    ctr_x= mean(mdl.nodes(:,1));
    ctr_y= mean(mdl.nodes(:,2));

% scale away from model

for e=1:length(mdl.electrode)
    if isfield(mdl.electrode(e),'nodes')
        elec_nodes= mdl.electrode(e).nodes;
        
        S= 1.00;
        vx= (mdl.nodes(elec_nodes,1) - ctr_x)*S;
        vy= (-mdl.nodes(elec_nodes,2) - ctr_y)*S; %still something wrong here
        % sort nodes around the model (to avoid crossed lines)
        [jnk,idx] = order_loop( [vx,vy], 0 );
    elseif isfield(mdl.electrode(e),'pos')
        vx = mdl.electrode(e).pos(:,1) - ctr_x;
        vy = mdl.electrode(e).pos(:,2) - ctr_y;
        idx = 1:length(vx);
    else
       eidors_msg('show_fem: WARNING: electrode %d has no nodes. Not showing.',e,2);
       continue;
    end
        
    ecolour = electr_colour( e );
    if numel(vx) == 1
       % Point Electrode Models: put a circle around the node
       line(vx(idx)+ctr_x,vy(idx)+ctr_y,  ...
            'LineWidth', 2, 'LineStyle','-','Color', ecolour, ...
            'Marker','o','MarkerSize', 6,'MarkerEdgeColor',ecolour);
    else
       % Complete/Shunt Electrode Models (multiple nodes per electrode)
       %  put a line along the edges that form the electrode
       line(vx(idx)+ctr_x,vy(idx)+ctr_y,  ...
            'LineWidth', 3, 'LineStyle','-','Color', ecolour, ...
            'Marker','none','MarkerSize', 6,'MarkerEdgeColor',ecolour);
    end
    if number_electrodes
       S= 1.05;
       vx= vx*S;
       vy= vy*S;
       switch number_electrodes
          case {1 true}
             txt = num2str(e);
          case 2
             try, txt = mdl.electrode(e).label; end
       end
       hh= text(mean(vx)+ctr_x, mean(vy)+ctr_y, txt);
       set(hh, 'HorizontalAlignment','center', 'FontWeight','bold');
    end
end
end

function ee= get_boundary( mdl )
   if isfield(mdl,'boundary')
       ee= mdl.boundary;
   else
       % calc and cache boundary
       ee = find_boundary( mdl.elems );
   end
end

function colour= electr_colour( e);
    if e==1;
       colour = [0,.7,0]; % light green electrode #1
    elseif e==2
       colour = [0,.5,0]; % mid-green electrode #2
    else
       colour = [0,.3,0]; % dark green
    end
end