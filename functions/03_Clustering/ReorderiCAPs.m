%% This function reorders the iCAP clusters according to the number of frames
% clusters
%
% IDX is the index of the cluster to which each trace belongs
%
% IDX2 is the rearranged index vector
function [IDX2,iCAPs2,iCAPs_folds2] = ReorderiCAPs(IDX,iCAPs,iCAPs_folds)

    IDX2 = [];
    n_clusters = max(IDX);
    
    % Computes how many traces of each type there are
    for iC = 1:n_clusters
        n_frames(iC) = sum(IDX==iC);
    end
    
    % Determines the order for which we have descending order of trace
    % amount
    [~,isort] = sort(n_frames,'descend');
    
    % Converts the idx vector so that the labels are now appropriate
    IDX2 = IDX;
    iCAPs2 = iCAPs;
    
    for iC = 1:n_clusters
        IDX2(IDX==isort(iC)) = iC;
        iCAPs2(iC,:) = iCAPs(isort(iC),:);
    end
    
    if exist('iCAPs_folds','var') && isstruct(iCAPs_folds)
        for iFold=1:length(iCAPs_folds.iCAPs) % reorder iCAPs of every fold
            % Now modifies the index and data appropriately
            for iC = 1:n_clusters
                IDX_new(iCAPs_folds.IDX{iFold}==isort(iC),1)=iC;
                iCAPs_new(iC,:)=iCAPs_folds.iCAPs{iFold}(isort(iC),:);
                iCAPs_z_new(iC,:)=iCAPs_folds.iCAPs_z{iFold}(isort(iC),:);
                sum_dist_new(iC,1)=iCAPs_folds.sum_dist{iFold}(isort(iC));
                dist_to_centroid_new(:,iC)=iCAPs_folds.dist_to_centroid{iFold}(:,isort(iC));
            end
            % save data on iCAPs folds
            iCAPs_folds2.iCAPs{iFold}=iCAPs_new;
            iCAPs_folds2.IDX{iFold}=IDX_new;
            iCAPs_folds2.sum_dist{iFold}=sum_dist_new;
            iCAPs_folds2.dist_to_centroid{iFold}=dist_to_centroid_new;
            iCAPs_folds2.iCAPs_z{iFold}=iCAPs_z_new;
        end
    else
        iCAPs_folds2=[];
    end
    
end