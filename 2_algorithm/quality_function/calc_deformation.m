
function deformation=calc_deformation(detected_organ,ref_organ)
% Function to calculate the deformation of the reconstructed organ

%ss62
    detected_organ=unique(detected_organ);
    length_intersection=length(intersect(detected_organ, ref_organ));
    deformation=(length(detected_organ)+length(ref_organ)-2*length_intersection)/length(ref_organ);
    %Should the deformation be truncated to 1 if it is greater than 1 ?
end

