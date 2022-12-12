%% This function performs the clustering of selected frames into a set of iCAPs
%
% Inputs:
% - Frames contains the matrix of frames on which we want to perform
% clustering (n_frames x n_ret_voxels)
% - K is the number of iCAPs into which we want to separate the data
% - DistType is a string defining the type of distance to use ('sqeuclidean'
% or 'cosine')
% - n_folds is a scalar defining how many times to run the k-means
% clustering process
function [iCAPs,IDX,dist_to_centroid,iCAPs_folds] = MakeiCAPs(Frames,param,fid)
    if ~exist(param.outDir_iCAPs,'dir');mkdir(param.outDir_iCAPs);end

    % in case the cluster distance to centroids should be saved for
    % evaluation, I need to write the loop over the several replicated
    % myself and keep all the results. Then the clusters will be matched to
    % each other using Hungarian algorithm (see Karahanoglu & Van De Ville,
    % EUSIPCO 2016) and the distances to the centroids across folds will be
    % saved in a cell array of length of the cluster number
    if ~isfield(param,'saveClusterReplicateData') || ~param.saveClusterReplicateData
        iCAPs_folds=[];
        if isfield(param,'MaxIter')
            [IDX,iCAPs,~,dist_to_centroid] = kmeans(Frames,param.K,'Distance',param.DistType,...
                'Replicates',param.n_folds,'Display','final',...
                'MaxIter',param.MaxIter);
        else
            [IDX,iCAPs,~,dist_to_centroid] = kmeans(Frames,param.K,'Distance',param.DistType,...
                'Replicates',param.n_folds,'Display','final');
        end
    else
        %% running clustering
        for iFold= 1:param.n_folds
            WriteInformation(fid,num2str(iFold))
            if isfield(param,'MaxIter')
                [IDX{iFold},iCAPs{iFold},sum_dist{iFold},dist_to_centroid{iFold}] = ...
                    kmeans(Frames,param.K,'Distance',param.DistType,...
                    'Replicates',1,'Display','final','MaxIter',param.MaxIter);
            else
                [IDX{iFold},iCAPs{iFold},sum_dist{iFold},dist_to_centroid{iFold}] = ...
                    kmeans(Frames,param.K,'Distance',param.DistType,...
                    'Replicates',1,'Display','final');
            end
        end
        
        %% match all clusters to results of first fold (Find corresponding
        % cluster centroids for every fold)
        disp('match all clustering results to first fold (Hungarian algorithm) ...')
        for iFold= 2:param.n_folds
            fprintf([num2str(iFold) ' '])
            dist_between_folds=pdist2(iCAPs{1},iCAPs{iFold},param.DistType);
            % hungarian
            [indexhun,costhun] = munkres(dist_between_folds);
            
            % reorder clusters
            for iC = 1:param.K
                IDX_new(IDX{iFold}==indexhun(iC),1)=iC;
                iCAPs_new(iC,:)=iCAPs{iFold}(indexhun(iC),:);
                sum_dist_new(iC,1)=sum_dist{iFold}(indexhun(iC));
                dist_to_centroid_new(:,iC)=dist_to_centroid{iFold}(:,indexhun(iC));
            end
            IDX{iFold}=IDX_new;
            iCAPs{iFold}=iCAPs_new;
            sum_dist{iFold}=sum_dist_new;
            dist_to_centroid{iFold}=dist_to_centroid_new;
            
            clearvars IDX_new iCAPs_new sum_dist_new dist_to_centroid_new
        end
        %% save data of iCAPs folds
        iCAPs_folds.iCAPs=iCAPs;
        iCAPs_folds.IDX=IDX;
        iCAPs_folds.sum_dist=sum_dist;
        iCAPs_folds.dist_to_centroid=dist_to_centroid;
        
        %% select best cluster (minimum total sum of distances)
        disp('Selecting best fold ...')
        for iFold= 1:param.n_folds
            fprintf([num2str(iFold) ' '])
            iCAPs_folds.total_dist_sum(iFold,1)=sum(iCAPs_folds.sum_dist{iFold});
            % z-scoring average iCAPs maps
            iCAPs_folds.iCAPs_z{iFold} = ZScore_iCAPs(iCAPs_folds.iCAPs{iFold},[],iCAPs_folds.IDX{iFold});
        end
        [~,bestID]=min(iCAPs_folds.total_dist_sum)
        iCAPs=iCAPs{bestID};
        IDX=IDX{bestID};
        sum_dist=sum_dist{bestID};
        dist_to_centroid=dist_to_centroid{bestID};
    end
    
    WriteInformation(fid,['iCAPs computed for ',num2str(param.K),' clusters',...
            ', ',num2str(param.n_folds),' folds and with distance ',param.DistType,'...']);
        
    
end