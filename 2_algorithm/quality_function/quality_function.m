function [quality_index_shape,quality_index_I,detailed_criteria,detailed_Dcond,coef_w,coef_bg] = quality_function(mode,img_ref,img_solve,options)
% mode = 1 calculates the quality index for absolute values of the
% reconstruction
% mode = 2 calculated the quality index for phase values of the
% reconstruction

%% Set variables
% org1 is lung, org2 is heart, org3 is spine, (org4 is fluid for case3)  
calc_modulus = options.calc_modulus; %ss2
relevant_cond= options.relevant_cond;%ss4
modified_norm= options.modified_norm;%ss6
npoints      = options.npoints; %ss8
maxit = options.maxit;
increment=0.025; options.increment=0.025;%ss15 % Increment for the window used in the organ detection (in percentage of the whole conductivity change range)
rec_call=false;%ss17 % Variable to detect whether the function was called recursively
disp_organs=options.disp_organs; %ss18% % to display the reference and detected organs at each iteration
% figure,fig=get(groot,'CurrentFigure'); idx_fig=fig.Number; %ss20 %ss22  Initialization of the number of the current figure
options.coef_quality=[]; coef_quality=[];% the coefficients for the computation of the overall quality   
bg_cond = options.relevant_cond_bg;


%% Normalise the variables
if mode ==1
  if calc_modulus
   relevant_cond=abs(relevant_cond); %ss28
   bg_cond=abs(bg_cond);
   img_ref.elem_data=abs(img_ref.elem_data);
   % [ref]cond-min(cond)
   relevant_cond=relevant_cond-min(img_ref.elem_data);
   bg_cond=bg_cond-min(img_ref.elem_data);

   img_ref.elem_data=img_ref.elem_data-min(img_ref.elem_data); %ss30
   
   % [ref]cond-min(cond)/(max(cond)-min(cond))
   relevant_cond=relevant_cond/max(img_ref.elem_data);
   bg_cond=bg_cond/max(img_ref.elem_data);
   img_ref.elem_data=img_ref.elem_data/max(img_ref.elem_data); %ss32

   % [solved]cond-min(cond)/(max(cond)-min(cond))
   if ~modified_norm %ss33
            img_solve.elem_data=abs(img_solve.elem_data);
            img_solve.elem_data=img_solve.elem_data-min(img_solve.elem_data);
            img_solve.elem_data=img_solve.elem_data/max(img_solve.elem_data);
    
   end
  end
end
if mode ==2
  if calc_modulus
   relevant_cond=radtodeg(angle((relevant_cond))); %ss28
   bg_cond=radtodeg(angle((bg_cond)));
   img_ref.elem_data=radtodeg(angle((img_ref.elem_data)));
   
   abs_ref = img_ref.elem_data;
   
   % [ref]cond-min(cond)
   relevant_cond=relevant_cond-min(img_ref.elem_data);
   bg_cond=bg_cond-min(img_ref.elem_data);
   img_ref.elem_data=img_ref.elem_data-min(img_ref.elem_data); %ss30
   
   % [ref]cond-min(cond)/(max(cond)-min(cond))
   relevant_cond=relevant_cond/max(img_ref.elem_data);
   bg_cond=bg_cond/max(img_ref.elem_data);
   img_ref.elem_data=img_ref.elem_data/max(img_ref.elem_data); %ss32

   % [solved]cond-min(cond)/(max(cond)-min(cond))
   if ~modified_norm %ss33
       img_solve.elem_data=radtodeg(angle((img_solve.elem_data)));
       
       % [solved] manage data
        phase_solve = img_solve.elem_data;
        elem_data_phase = phase_solve-min(phase_solve);
        elem_data_phase = elem_data_phase/max(elem_data_phase);

        elem_data_phase = elem_data_phase.*(max(abs(abs_ref)));
        elem_data_phase = elem_data_phase+min(abs_ref);
        img_solve.elem_data = elem_data_phase;
       
       img_solve.elem_data=img_solve.elem_data-min(img_solve.elem_data);
       img_solve.elem_data=img_solve.elem_data/max(img_solve.elem_data);
   end
  end
end  
%% Calculation of the radius of the domain %%ss38
% 0.07,1.6243,1.38
radius_min=min(sqrt(img_ref.fwd_model.nodes(img_ref.fwd_model.boundary(:,2),1).^2+img_ref.fwd_model.nodes(img_ref.fwd_model.boundary(:,2),2).^2));
radius_max=max(sqrt(img_ref.fwd_model.nodes(img_ref.fwd_model.boundary(:,2),1).^2+img_ref.fwd_model.nodes(img_ref.fwd_model.boundary(:,2),2).^2));
radius_pix=radius_min/radius_max*npoints/2;

%% Create pixel grid
    % ss39
img_ref.calc_colours.npoints=npoints;
ref_grid=calc_slices(img_ref);
size_domain=npoints*npoints-sum(sum(isnan(ref_grid)));
img_ref=rmfield(img_ref,'calc_colours'); %This line to ensure that the organs are well displayed

img_solve.calc_colours.npoints=npoints;
recons_grid=calc_slices(img_solve);
    
figure,subplot(1,2,1),imshow(ref_grid),title('ref\_grid'),
subplot(1,2,2),imshow(recons_grid),title('recons\_grid'),

%% Organ detection and computation of the quality index
quality_index_shape=0; quality_index_I=0;
detailed_quality=zeros(length(relevant_cond),1);
detailed_Dcond = zeros(length(relevant_cond),1);
detailed_criteria=zeros(4,length(relevant_cond));
sum_coef=0; coef_w=[];
idx_fig=200;
for idx1=1:length(relevant_cond) %for loop where the quality index for each organ is calculated
    if idx1 ==1
        max_iter = maxit(1);
    end
    if idx1 ==2
        max_iter = maxit(2);        
    end
    if idx1 ==3
        max_iter = maxit(3);        
    end
    if idx1 ==4
        max_iter = maxit(4);       
    end
    % 3 organs    
        % ss42
        quality_organ=Inf; %Inf is the baddest possible quality
        Dcond_organ = Inf; %Inf is the baddest possible quality
        
        % The black and white image of the organ is created. Then the
        % number of organs with the given conductivity is calculated (important for the case of the lungs)
        organ_grid=abs(ref_grid-relevant_cond(idx1))<=0.01;%eps;
        
        %Case where no organ with the given conductivity was detected (i.e. a false value was given in relevant_cond)

%         if sum(sum(organ_grid))<eps
%             %ss43xx
%             detailed_quality(idx1)=NaN;
%             continue
%         end
        
        if isempty(coef_quality)
            coef=sum(sum(organ_grid));
            %ss44
        else
            coef=coef_quality(idx1);
            %ss45xx
        end
        sum_coef=sum_coef+coef;
         %ss46
        if disp_organs
             %ss47xx
            figure(idx_fig+1);
            img_test_ref = calc_colours(organ_grid,img_ref);
            img_test_ref(organ_grid==0)=128;
            img_test_ref(organ_grid==1)=255;
            subplot(1,length(relevant_cond),idx1);
            image(img_test_ref);
            axis square
        end
        %ss48
        detection_infos_ref=bwconncomp(organ_grid,8);
        num_organs=detection_infos_ref.NumObjects;
        %ss59,ss60
        idx_organs_sorted=find_max_cell(detection_infos_ref.PixelIdxList,num_organs); 
        %indices of the organs in the reference in the list PixelIdxList, 
        %sorted by associating each reference organ with the nearest detected one
        
        center_ref=zeros(num_organs,2);
        for idx3=1:num_organs %Calculation of the center of the organs in the reference image
            %ss49
            %ss58
           center_ref(idx3,:)=calc_center(detection_infos_ref.PixelIdxList{idx_organs_sorted(idx3)},npoints);
        end                    

        %--------------------------For loop to calculate the quality index for each organ/relevant conductivity. 
        % Only the best quality is conserved
        for idx2=1:max_iter  %each adaptation of ref-window on recons_grid
             %ss50
            min_wind=relevant_cond(idx1)-idx2*increment; %min cond of window
            max_wind=relevant_cond(idx1)+idx2*increment; %max cond of window
            %ss56
            [organ_bw, detection_infos]=organ_detection(recons_grid,min_wind,max_wind,organ_grid);
            
           
            if disp_organs
                 %ss51xx
                figure(idx_fig+1+idx1);
                img_test_ref = calc_colours( organ_bw, img_ref);
                img_test_ref(organ_bw==0)=128;
                img_test_ref(organ_bw==1)=255; 
                subplot(1,max_iter,idx2);
                image(img_test_ref);
                axis square
            end
            
            if detection_infos.NumObjects>0 %We only procedd if an organ was detected otherwise, it is not necessary to update the quality_organ and we can proceed to the next iteration 
                 %ss52
                 %ss64,ss58,ss57,ss65
                idx_detected_organs=find_organs_ref(detection_infos.PixelIdxList,center_ref,npoints);
                %Calculation of the displacement and deformation for the
                %num_organs nearest from the reference detected organs

                total_shape_ref=[]; % For the computation of the shape deformation, we consider the organs as a whole
                for idx=1:num_organs
                     %ss53
                     %sorted vectors of cond for ref
                    total_shape_ref=[total_shape_ref;detection_infos_ref.PixelIdxList{idx_organs_sorted(idx3)}];
                end
                
                total_displacement=0;
                total_shape_organs=[]; % For the computation of the shape deformation, we consider the detected organs as whole
                num_iter=length(idx_detected_organs); % The number of iterations is the number of detected organs in the reconstruction
                for idx3=1:num_iter
                     %ss54
                     %ss58,57
                    calc_center(detection_infos.PixelIdxList{idx_detected_organs(idx3)},npoints);
                    %ss61,58,57
                    total_displacement=total_displacement+...
                        calc_displacement(detection_infos.PixelIdxList{idx_detected_organs(idx3)},center_ref(idx3,:),npoints,radius_pix);
                    total_shape_organs=[total_shape_organs;detection_infos.PixelIdxList{idx_detected_organs(idx3)}];
                end
                    
                    total_displacement=(total_displacement+(num_organs-num_iter))/num_organs; % We add a penalty if fewer organs were detected than there are in the reference
                    %ss62
                    total_shape_def=calc_deformation(total_shape_organs,total_shape_ref)/num_organs; % We calculate the shape deformation by considering the organs in the ref and reconstruction as a whole
                    %ss63
                    penalty=calc_penalty(detection_infos.PixelIdxList,unique(idx_detected_organs),total_shape_ref); 
                    
                    
                    quality_organ_iter=(total_displacement+total_shape_def)/2+penalty;               

                if quality_organ>quality_organ_iter
                     %%ss55
                    quality_organ=quality_organ_iter;
                    detailed_criteria(:,idx1)=[total_displacement;total_shape_def;penalty;idx2];
                end 
            end %end of if detected organ
            
        end
        
        
        % intersection btw recons, ref
            selected_bw = bitand(organ_bw,organ_grid);
            
            selected_indx{idx1}.recons = find(organ_bw==1);
            selected_indx{idx1}.ref = find(organ_grid==1);
            selected_indx{idx1}.intesect = find(selected_bw==1);
           % assign normalized cond to intersected region
           selected_recons_org = zeros(64,64);
            selected_recons_org(selected_indx{idx1}.intesect) = recons_grid(selected_indx{idx1}.intesect);
        
         % cond diff for each organ
        selected_ref_org = zeros(npoints,npoints);
        if isempty(selected_indx{idx1}.intesect) || detailed_criteria(3,idx1) ~= 0
            selected_ref_org(selected_indx{idx1}.ref) = 1;
        else
            selected_ref_org(selected_indx{idx1}.intesect) = relevant_cond(idx1);
        end
        selected_ref_org_disp = zeros(npoints,npoints);
        selected_ref_org_disp(selected_indx{idx1}.ref) = relevant_cond(idx1);
            
        % display selected organs
        if disp_organs
            figure,subplot(3,2,1),imshow(organ_bw),title('recons organ indx1'),
                   subplot(3,2,2),imshow(organ_grid),title('ref organ indx1');
                   subplot(3,2,3),imshow(selected_bw),title('intersected organ(bw) indx1');
            subplot(3,2,4),imshow(selected_recons_org),title('intersected organ(grey) indx1');
            subplot(3,2,5),imshow(selected_ref_org);title('ref cond, intersected'),
            subplot(3,2,6),imshow(selected_ref_org_disp);title('ref cond, intersected disp'),
        end
             
        if detailed_criteria(3,idx1) ~= 0 %if no detected organ (Penalty ~=0)
            Dcond_organ = 1;
        else
            Dcond_organ = calc_cond_diff(selected_ref_org, selected_recons_org,coef);
        end
        detailed_Dcond(idx1) = Dcond_organ;
        
        % calculate detailed_Q,Q-index 
        detailed_quality(idx1)=quality_organ;
        
        % calculate final shape based, intensity based(organs) Q_index
        quality_index_shape = quality_index_shape + (coef*quality_organ);
        quality_index_I = quality_index_I + (coef*Dcond_organ);
        coef_w(idx1) = coef;

end

%% bg conductivity
ref_all_org_idx = []; recons_all_org_idx = [];
for i=1:idx1
 ref_all_org_idx = [ref_all_org_idx; selected_indx{i}.ref];
 recons_all_org_idx = [recons_all_org_idx; selected_indx{i}.intesect];
end
ref_all_org_idx = unique(ref_all_org_idx); recons_all_org_idx = unique(recons_all_org_idx);
ref_bg_idx = (setdiff([1:npoints*npoints], ref_all_org_idx));
ref_bg_grid = zeros(npoints,npoints); ref_bg_grid(ref_bg_idx)=ref_grid(ref_bg_idx);
recons_bg_idx = setdiff(setdiff([1:npoints*npoints], recons_all_org_idx),ref_all_org_idx);
recons_bg_grid = zeros(npoints,npoints); recons_bg_grid(recons_bg_idx)=recons_grid(recons_bg_idx);

ref_bg_grid(isnan(recons_bg_grid) == 1) = NaN; %1 is NaN % manage bg area of reference 
recons_bg_grid(isnan(ref_bg_grid) == 1) = NaN; %1 is NaN % manage bg area of reconstruction 

coef_bg = length(find(isnan(ref_grid(ref_bg_idx))==0));
if disp_organs
    figure,subplot(1,2,1),imshow(ref_bg_grid),title('conductivity on ref bg');
    subplot(1,2,2),imshow(recons_bg_grid),title('conductivity on recons bg');
end
Dcond_bg = calc_cond_diff(ref_bg_grid,recons_bg_grid,coef_bg);
detailed_Dcond(idx1+1) = Dcond_bg;
quality_index_I = (quality_index_I + coef_bg*Dcond_bg);

quality_index_I = quality_index_I/(sum(coef_w)+coef_bg);
quality_index_shape = quality_index_shape/(sum(coef_w));