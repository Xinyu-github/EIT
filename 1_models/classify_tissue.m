function [pixelmap_img_class, area_indexes, safety_indexes_map, safety_indexes, classification_map] = classify_tissue(input_img, freq_1, freq_2, res, contour_levels_phase, contour_levels_mag, use_mag_or_phase, show_boarders, show_steps)
% This function classifies different tissues based on their 
% differential conductivities as taken from an EIDORS image. The function
% is ideally used on images that were reconstructed using the
% mixed-minimization algorithm developed at MedIT.
% The algorithm can currently differentiate between trabecular bone, 
% cortical bone, blood, heart, inflated lung, deflated lung, muscle and
% fluid. This function is NOT part of the EIDORS package and was developed 
% at the MEDIT institute.
%
% Syntax: [pixelmap_img_class, area_indexes, safety_indexes_map, safety_indexes, classification_map] = classify_tissue(input_img, freq_1, freq_2, res, contour_levels_phase, contour_levels_mag, use_mag_or_phase, show_boarders, show_steps)
%
% Input arguments:
%   input_img:              An EIDORS image that represents an dfEIT image. 
%                           Important: The image object has to consist of 
%                           differential conductivities, otherwise the 
%                           classification will yield wrong results.
%   freq_1, freq_2:         The two frequencies at which the measurements 
%                           for the dfEIT image were performed.
%   phase_round:            Factor by which the phase values of the dfEIT 
%                           image should be rounded. The rounding factir 
%                           defines the sensitivity with which organs are 
%                           separated. 
%   res:                    Resolution of the final classified image.
%   contour_levels_phase:   Defines how many contour lines will be
%                           extracted from the phase data.
%   contour_levels_mag:     Defines how many contour lines will be
%                           extracted from the magnitude data.
%   use_mag_or_phase:       Defines whether only phase (1), only magnitude
%                           (2) or the fusion of both information (3) shall
%                           be used for segmentation.
%   show_boarders:          Boolean that defines whether the segmentation
%                           lines are present in the classified image. If
%                           true, boarder-pixels will have the value '0'.
%   show_steps:             Boolean that defines whether the different 
%                           steps of image processing shall be plotted. 
%
% Output arguments:
%   pixelmap_img_class: A matrix of size resXres. This matrix displays the
%                       classified EIT image in pixel form rather than as
%                       an EIDORS image object. Every pixel is given a
%                       complex conductivity value depending on the tissue
%                       type it was classified as.
%   area_indexes:       A matrix of size resXres. The matrix represents the
%                       result of the segmentation process. Every segment
%                       is defined by a indivudutal integer starting at 1
%                       and then ascending up to the number of total
%                       segments. The boarders between the segments are
%                       makred by the number 0. The image can be displayed
%                       as an image if desired.
%   safety_indexes_map: A matrix of size resXres that specifies the safety
%                       factor for each segment. Each pixel in each segment
%                       is given the safety factor of the segment in %. The
%                       boarders between the segments are given the number
%                       0.
%   safety_indexes:     An array of all the safety indexes for the
%                       different segments. The index of each factor
%                       corresponds to the number of the given segment as
%                       defined in the area_indexes matrix.
%   classification_map: A matrix of size resXres. This matrix displays the
%                       classified EIT image in pixel form rather than as
%                       an EIDORS image object. Every pixel is given a
%                       number that specifies the tissue type it was
%                       classified as. Tissue types are specified with 
%                       numbers from 1 to 8 as follows: blood (1), 
%                       cancellous bone(2), cortical bone (3), heart (4),
%                       inflated lung (5), deflated lung (6),
%                       muscle (7), fluid (8)



% Initialise the classified image and the array that stores organ positions
img_class = input_img;

%% Load the conductivity reference data

% Load the tissue conductivity data
tissues_f1 = getTissues_mod(freq_1);
tissues_f2 = getTissues_mod(freq_2);

% Calculate the frequency-differential conductivities
tiss_diff.blood = round(tissues_f2.Blood_combined - tissues_f1.Blood_combined, 4, 'decimal');
tiss_diff.bone_schwamm = round(tissues_f2.Bone_Schwamm_combined - tissues_f1.Bone_Schwamm_combined, 4, 'decimal');
tiss_diff.bone_outside = round(tissues_f2.BoneAussen_combined - tissues_f1.BoneAussen_combined, 4, 'decimal');
tiss_diff.heart = round(tissues_f2.Heart_combined - tissues_f1.Heart_combined, 4, 'decimal');
tiss_diff.lung_in = round(tissues_f2.LungInflated_combined - tissues_f1.LungInflated_combined, 4, 'decimal');
tiss_diff.lung_de = round(tissues_f2.LungDeflated_combined - tissues_f1.LungDeflated_combined, 4, 'decimal');
tiss_diff.muscle = round(tissues_f2.Muscle_combined - tissues_f1.Muscle_combined, 4, 'decimal');
tiss_diff.fluid = round(tissues_f2.fluid_combined - tissues_f1.fluid_combined, 4, 'decimal');

% calculate the magnitude of the stored conductivities
tiss_abs_diff = structfun(@abs, tiss_diff);

%% Create separated areas in the image to distinguish between organs and apply the best fitting tissue type to every pixel

% Use only magnitude data for area segmentation
if use_mag_or_phase == 1
    
    % Run separate_phase_areas to get the different areas that are separated by
    % phase changes
    [area_indexes, area_number] = separate_magnitude_areas(img_class, res, contour_levels_mag, show_steps);
    
    if show_steps
        figure('Name', 'Segmentation Map', 'Position', [1150 50 500 400]);
        imshow(area_indexes,[],'InitialMagnification','fit');
    end
    
    % Get the magnitude of the input image
    img_class.elem_data = abs(img_class.elem_data);
    pixelmap_img_class = calc_slices(img_class);
    delta_cond = zeros(length(tiss_abs_diff), 1);
    
% Use only phase data for area segmentation
elseif use_mag_or_phase == 2
        
    % Run separate_phase_areas to get the different areas that are separated by
    % phase changes
    [area_indexes, area_number] = separate_phase_areas(img_class, res, contour_levels_phase, show_steps);
    
    if show_steps
        figure('Name', 'Segmentation Map', 'Position', [1150 50 500 400]);
        imshow(area_indexes,[],'InitialMagnification','fit');
    end
    
    % Get the magnitude of the input image
    img_class.elem_data = abs(img_class.elem_data);
    pixelmap_img_class = calc_slices(img_class);
    delta_cond = zeros(length(tiss_abs_diff), 1);
     
% Use fused information from magnitude and phase for area segmentation
elseif use_mag_or_phase == 3

    % Run separate_phase_areas to get the different areas that are separated by
    % phase changes
    [area_indexes_phase, area_no_phase] = separate_phase_areas(img_class, res, contour_levels_phase, show_steps);
    [area_indexes_mag, area_no_mag] = separate_magnitude_areas(img_class, res, contour_levels_mag, show_steps);
 
    area_indexes = ones(size(area_indexes_phase));
    area_indexes(area_indexes_phase == 0) = 0;
    area_indexes(area_indexes_mag == 0) = 0;
    area_indexes = bwlabel(area_indexes);
    
    if show_steps
        area_indexes_show = area_indexes;
        area_indexes_show(area_indexes_show ~= 0) = 1;
        
        figure('Name', 'Segmentation Map', 'Position', [1150 50 500 400]);
        imshow(area_indexes_show,[],'InitialMagnification','fit');
    end
    
    % Get the magnitude of the input image
    img_class.elem_data = abs(img_class.elem_data);
    pixelmap_img_class = calc_slices(img_class);
    tissue_no = length(tiss_abs_diff);
    delta_cond = zeros(tissue_no, 1);
    area_number = max(max(area_indexes));
       
end

% set boards to zero
pixelmap_img_class((area_indexes == 0) & ~isnan(pixelmap_img_class)) = 0;

% define matrix to store the safety indexes and the classified image
safety_indexes_map = zeros(res, res);
safety_indexes = zeros(length(area_number));
classification_map = zeros(res, res);

% Iterate over all phase areas and find the most fitting tissue type for
% each area based on the average magnitude within each area.
for k = 1:area_number
       
    % Determine the total number of pixels in this segment
    pixels_in_segment = sum(sum(area_indexes == k));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Firstly: classify entire area as one tissue type
    mean_abs = mean(pixelmap_img_class(area_indexes == k));
    
    % Find the delta between the current average magnitude and the
    % magnitude of each tissue type.
    for j = 1:8
        delta_cond(j) = mean_abs - tiss_abs_diff(j);
    end
    
    % Determine the minimum delta and thus the best fitting tissue type for
    % the current phase area.
    [min_delta, idx_min_delta] = min(abs(delta_cond));
    classification_map(area_indexes == k) = idx_min_delta;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Secondly: Calculate Safety-Factor for each segment
        correct_pixel_counter = 0; 
        for i = 1:res
            for j = 1:res
                % Repeat the classification procedure for every pixel in
                % the current segment
                if area_indexes(i,j) == k          
                    delta_cond = zeros(tissue_no, 1);
                    current_pixel_abs = pixelmap_img_class(i,j);
                    for l = 1:8
                        delta_cond(l) = current_pixel_abs - tiss_abs_diff(l);
                    end
                    [min_delta, current_pixel_idx_min_delta] = min(abs(delta_cond));
                    % Count the number of pixels whose classification
                    % coincides with the classification of the segment
                    if current_pixel_idx_min_delta == idx_min_delta
                        correct_pixel_counter = correct_pixel_counter + 1;
                    end 
                end
            end
        end
        % Calculate the safety factor of the current segment and store it
        % in the safety_indexes pixelmap
        safety_factor = (correct_pixel_counter/pixels_in_segment)*100;
        safety_indexes_map(area_indexes == k) = safety_factor;    
        safety_indexes(k) = safety_factor;
end

% Convert area indexes to double in case it wasn't already to get a
% consistend output parameter
area_indexes = double(area_indexes);

% If show_boarders is false, remove the black boarderlines by replacing
% them with the most common classification in the surrounding pixels that
% is not zero.
if ~show_boarders
    for i = 1:res
        for j = 1:res
            if pixelmap_img_class(i,j) == 0
                adjacent = classification_map(i-1:i+1, j-1:j+1);
                replacement = mode(mode(adjacent));
                if replacement == 0
                    adjacent(adjacent == 0) = NaN;
                    classification_map(i,j) = mode(mode(adjacent));
                else
                    classification_map(i,j) = replacement;
                end
            end
        end
    end    
end

% Plot the classification image
plot_classification_map(classification_map, safety_indexes_map);

%% Map the complex conductivities

pixelmap_img_class(classification_map == 1) = tiss_diff.blood;
pixelmap_img_class(classification_map == 2) = tiss_diff.bone_schwamm;
pixelmap_img_class(classification_map == 3) = tiss_diff.bone_outside;
pixelmap_img_class(classification_map == 4) = tiss_diff.heart;
pixelmap_img_class(classification_map == 5) = tiss_diff.lung_in;
pixelmap_img_class(classification_map == 6) = tiss_diff.lung_de;
pixelmap_img_class(classification_map == 7) = tiss_diff.muscle;
pixelmap_img_class(classification_map == 8) = tiss_diff.fluid;

end
