function idx_organs_sorted=find_organs_ref(detected_organs,center_ref,npoints)
% Function to associate each reference with its nearest detected organ.
    
%ss64
    aux_fun=@(x) calc_center(x,npoints);
    centers=reshape(cell2mat(cellfun(aux_fun,detected_organs,'UniformOutput',false)),2,[])';
    centers=centers(:,1)+1i*centers(:,2);
    center_ref=center_ref(:,1)+1i*center_ref(:,2);

    idx_organs_sorted=zeros(length(center_ref),1);
    for idx=1:length(center_ref)   
        %ss65
       [~,idx_opt]=min(abs(centers-center_ref(idx)));
       % One detected organ can be taken twice. It is not a problem
       idx_organs_sorted(idx)=idx_opt;
    end


end

