function penalty=calc_penalty(detected_organs_list, idx_detected_organs,total_shape_ref)
%Function to compute the penalty of the organ detection. This penalty is
%the number of the too many detected organs times the number of pixels that
%we detected but are not part of the organ over the total number of pixels
%in the image

%ss63
aux=cellfun(@length, detected_organs_list);
penalty=(length(aux)-length(unique(idx_detected_organs)))*(sum(aux)-sum(aux(unique(idx_detected_organs))))/length(total_shape_ref);

end

