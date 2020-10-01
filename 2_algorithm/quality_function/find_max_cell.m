function idx_max=find_max_cell(cell_array,num_max)
% ss59
%Function to find the first num_max biggest elements of a cell array sorted
%by growing sizes
if num_max==0
   %ss60
    idx_max=[];
    return
end

aux=cellfun(@length,cell_array);
[~,idx]=max(aux);
cell_array{idx}=[];
idx_max=[idx,find_max_cell(cell_array,num_max-1)];

end
