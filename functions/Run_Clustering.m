%% This function runs clustering on the significant innovations found 
%  previously by thresholding, results will be the iCAPs maps
%
% Input: 
%   param - struct containing all necessary parameters to run TA
% 
%     * Data reading and saving information:
%       .PathData - path to data
%       .Subjects - cell array with list of subdirectories where each 
%           subject's fMRI data is stored; must contain one entry per 
%           subject to analyze.
%           This is where the TA folder will be created (or looked for) 
%           for each subject.
%       .n_subjects - number of subject to analyze
%       .title - possibility to define a title for the current
%           project (usefull if TA should be run for different 
%           parameters), default: current date, it should be the title of
%           the project for which TA has already been run
%       [.force_Aggregating] - if set to 1, significant innovations
%           and activity-inducing signal will be read and aggregated, even
%           if this has been done before
%       [.force_Clustering] - if set to 1, Clustering will be forced to 
%           run, even if already has been done
%       [.thresh_title] - information to read thresholding data and save
%           clustering results, if not specified, the fields 'param.alpha'
%           and 'param.f_voxels' need to exist
%       [.data_title] - information for saving of aggregated clustering data
%       [.iCAPs_title] - string or cellstring with all the subfolders to
%           create for saving clustering results
% 
%     * Clustering information:
%       [.doClustering] - specify if clustering should be done (one may
%           want to only run consensus clustering, then set this to 0),
%           default = 1
%       [.saveClusterReplicateData] - specify if the result of each replicate 
%           should be saved during clustering, default = 0
%       .n_folds - Number of replicates of clustering
%       .K - Number of clusters, can be a number or an array of multiple K
%       .DistType - Type of distance to use for the k-means clustering 
%           process (choose between 'sqeuclidean' and 'cosine')
%       [.MaxIter] - Maximum number of allowed iterations of the kmeans 
%           clustering, the Matlab default of 100 is sometimes not enough 
%           if many frames (i.e., many subjects or long scans) are 
%           included, default = 100
% 
%     * Consensus clustering information:
%       [.doConsensusClustering] - specify if consensus clustering should
%           be done on top of the clustering as specified above
%       [.force_ConsensusClustering] - if set to 1, Consensus Clustering 
%           will be forced to run, even if already has been done
%       .Subsample_type - Subsampling Type of consensus clustering
%           (default: 'items')
%           'subjects' to subsample all frames from a subject; 
%           'items' to subsample frames without taking into account 
%               within- or between-subject information
%       .Subsample_fraction - fraction of subsampled data per fold
%       .cons_n_folds - number of folds for consensus clustering
%           (clustering will be run with 'param.n_folds' replicates in each
%           consensus clustering fold)
%       [.cons_title] - subfolder in which to save consensus results
%
% Output:
%   Creates iCAPs maps and consensus clustering results


function [] = Run_Clustering(param)

    % Date and time when the routines are called
    param.date = strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
    if ~isfield(param,'title') || isempty(param.title)
        param.title=param.date;
    end
    
    % Creates and opens a log-file that will contain all information related to
    % what is done to the data
    if ~exist(fullfile(param.PathData,'TAlogs'),'dir'); mkdir(fullfile(param.PathData,'TAlogs'));end;
    fid=fopen(fullfile(param.PathData,'TAlogs',['log_Clustering_',param.title,'.txt']),'a+');
    WriteInformation(fid,['Starting the total activation/iCAPs tools for project entitled ',param.title]);
    
    % Checks that the path towards the data is correct
    if ~exist(param.PathData,'dir')
        WriteInformation(fid,'Incorrect path towards the data: execution stopped');
        error('The data folder that you specified does not exist ! Please check and restart running...');
    end
    
    % setting data title, if not specified
    if ~isfield(param,'data_title') || isempty(param.data_title)
        param.data_title=param.title;
    end
    
    % setting thresholding title, if not specified
    if ~isfield(param,'thresh_title') || isempty(param.thresh_title)
        param.thresh_title = ['Alpha_',strrep(num2str(param.alpha(1)),'.','DOT'),'_',...
            strrep(num2str(param.alpha(2)),'.','DOT'),'_Fraction_',...
            strrep(num2str(param.f_voxels),'.','DOT')];
    end
    
    % main folder according to data+thresholding, aggregated frames will be
    % saved here
    param.outDir_main=fullfile(param.PathData,'iCAPs_results',[param.data_title,'_',param.thresh_title]);
    if ~exist(param.outDir_main,'dir');mkdir(param.outDir_main);end;
    
    % setting iCAPs title (subfolder in which to save clustering results)
    if ~isfield(param,'iCAPs_title') || isempty(param.iCAPs_title)
        for nK=1:length(param.K)
            param.iCAPs_title{nK} = ['K_',num2str(param.K(nK)),'_Dist_',...
                    param.DistType,'_Folds_',num2str(param.n_folds)];
        end
        if length(param.K)==1
            param.iCAPs_title=param.iCAPs_title{1};
        end
    end
    
    
    
    % check if clustering has already been done
    [Aggregating_done,~,~] = Check_iCAPs_Files(param.outDir_main);
    if isfield(param,'force_Aggregating') && param.force_Aggregating
        Aggregating_done=0;
    end
    
    % only load aggregated data if clustering needs to be done (if only consensus is to be calculated, we don't need AI etc)
    if (isfield(param,'doConsensusClustering') && param.doConsensusClustering) || (~isfield(param,'doClustering') || param.doClustering)
    if ~Aggregating_done
        % aggregating frames
        WriteInformation(fid,'Aggregating Subject Frames...');
        
        % TA and thresholding paths for all subjects
        SubjPath = fullfile(param.PathData,param.Subjects);
        resultsPath=fullfile(SubjPath,'TA_results',param.title);
        
        [I_sig,final_mask,subject_labels,time_labels,AI,AI_subject_labels] = AggregateSubjectFrames(resultsPath,param,fid);
        WriteInformation(fid,'Saving aggregated data...');
        save(fullfile(param.outDir_main,'I_sig.mat'),'I_sig','-v7.3');
        save(fullfile(param.outDir_main,'AI.mat'),'AI','-v7.3');
        save(fullfile(param.outDir_main,'AI_subject_labels.mat'),'AI_subject_labels','-v7.3');
        save(fullfile(param.outDir_main,'subject_labels.mat'),'subject_labels','-v7.3');
        save(fullfile(param.outDir_main,'time_labels.mat'),'time_labels','-v7.3');
        save(fullfile(param.outDir_main,'final_mask.mat'),'final_mask','-v7.3');
        save4Dnii(param.outDir_main,'','final_mask',final_mask,fullfile(resultsPath{1},'Thresholding',param.thresh_title,'mask_nonan.nii'));
    else
        WriteInformation(fid,'Aggregating already done, loading aggregated data...');
        % loading aggregated data
        load(fullfile(param.outDir_main,'I_sig.mat'));
        load(fullfile(param.outDir_main,'final_mask.mat'));
        load(fullfile(param.outDir_main,'subject_labels.mat'));
        load(fullfile(param.outDir_main,'time_labels.mat'));
        load(fullfile(param.outDir_main,'AI.mat'));
        load(fullfile(param.outDir_main,'AI_subject_labels.mat'));
    end
    end
    
    % adapting k and iCAPs title to make the code compatible with
    % multiple K
    if iscell(param.iCAPs_title)
        param.iCAPs_title_cell=param.iCAPs_title;
    else
        param.iCAPs_title_cell{1}=param.iCAPs_title;
    end
    param.K_vect=param.K;
    
    % main procedure for clustering:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % consensus clustering
    if isfield(param,'doConsensusClustering') && param.doConsensusClustering
        % setting name of subfolder in which to save consensus results
        if ~isfield(param,'cons_title') || isempty(param.cons_title)
            param.cons_title=[num2str(param.K(1)) 'to' num2str(param.K(end)) ...
                '_SubsampleType_' param.Subsample_type ...
                '_Fraction_' strrep(num2str(param.Subsample_fraction),'.','DOT') ...
                '_nFolds_' num2str(param.cons_n_folds) ...
                '_Dist_' param.DistType];
        end
        
        % defining consensus clustering output directory
        param.outDir_cons=fullfile(param.PathData,'iCAPs_results',...
            [param.data_title,'_',param.thresh_title],param.cons_title);
        if ~exist(param.outDir_cons,'dir');mkdir(param.outDir_cons);end;
        
        
        % check if consensus clustering has already been done
        [~,~,ConsensusClustering_done,~] = Check_iCAPs_Files([],[],param.outDir_cons);
        
        if isfield(param,'force_ConsensusClustering') && param.force_ConsensusClustering
            ConsensusClustering_done=0;
        end
        
        if ~ConsensusClustering_done
            ConsensusClustering(I_sig,subject_labels,param);
            save(fullfile(param.outDir_cons,'param'),'param','-v7.3');
        end
        
    end
    
    % clustering
    if ~isfield(param,'doClustering') || param.doClustering
        WriteInformation(fid,'Entering the clustering process...');
        % clustering for every K
        for iK=1:length(param.iCAPs_title_cell)

            param.iCAPs_title=param.iCAPs_title_cell{iK};
            param.K=param.K_vect(iK);

            WriteInformation(fid,['K = ' num2str(param.K) ', iCAPs title: ' param.iCAPs_title]);
            
            % clustering subfolder to save clustering-specific information
            param.outDir_iCAPs=fullfile(param.outDir_main,param.iCAPs_title);
            if ~exist(param.outDir_iCAPs,'dir');mkdir(param.outDir_iCAPs);end;
            
            % check if clustering has already been done
            [~,Clustering_done,~] = Check_iCAPs_Files(param.outDir_main,param.outDir_iCAPs);
            
            if isfield(param,'force_Clustering') && param.force_Clustering
                Clustering_done=0;
            end
            
            if ~Clustering_done
                % Clustering of the iCAPs
                WriteInformation(fid,'Running Clustering...');
                [iCAPs,IDX,dist_to_centroid,iCAPs_folds] = MakeiCAPs(I_sig,param,fid);

                % Rearranges the data and centroids and indices so that '1' denotes
                % the cluster with most occurrences in the data
                WriteInformation(fid,'Rearranging Time Courses...');
    %                 [Data,IDX,iCAPs,subject_labels,time_labels,iCAPs_folds] = RearrangeTimeCourses(Data,IDX,iCAPs,subject_labels,time_labels,iCAPs_folds);
                % Dani: I replaced this by a function which just re-orders
                % the iCAPs (i.e. iCAPs numbering according to innovation 
                % frame counts), but without changing the order in the Data
                [IDX,iCAPs,iCAPs_folds] = ReorderiCAPs(IDX,iCAPs,iCAPs_folds);

                %% z-score iCAPs
                WriteInformation(fid,'z-scoring iCAPs...');
                iCAPs_z = ZScore_iCAPs(iCAPs,I_sig,IDX);

                if exist('iCAPs_folds','var') && ~isempty(iCAPs_folds)
                    for iFold=1:length(iCAPs_folds.iCAPs)
                        iCAPs_folds.iCAPs_z{iFold} = ZScore_iCAPs(iCAPs_folds.iCAPs{iFold},I_sig,iCAPs_folds.IDX{iFold});
                    end
                end

                %% saving
                WriteInformation(fid,'Saving clustering data...');
                save(fullfile(param.outDir_iCAPs,'iCAPs'),'iCAPs','-v7.3');
                save(fullfile(param.outDir_iCAPs,'iCAPs_z'),'iCAPs_z','-v7.3');
                save(fullfile(param.outDir_iCAPs,'IDX'),'IDX','-v7.3');
                save(fullfile(param.outDir_iCAPs,'dist_to_centroid'),'dist_to_centroid','-v7.3');
                save(fullfile(param.outDir_iCAPs,'param'),'param','-v7.3');
                if ~isempty(iCAPs_folds)
                    save(fullfile(param.outDir_iCAPs,'iCAPs_folds'),'iCAPs_folds','-v7.3');
                elseif exist(fullfile(param.outDir_iCAPs,'iCAPs_folds.mat'),'file')
                    delete(fullfile(param.outDir_iCAPs,'iCAPs_folds.mat'));
                end
                save4Dnii(param.outDir_iCAPs,'','iCAPs',iCAPs',fullfile(param.outDir_main,'final_mask.nii'),final_mask);
                save4Dnii(param.outDir_iCAPs,'','iCAPs_z',iCAPs_z',fullfile(param.outDir_main,'final_mask.nii'),final_mask);
                
            else
                 WriteInformation(fid,'Clustering already done, skipping...');
            end
            
            %% saving subjects maps
            subjectsSaved=exist(fullfile(param.outDir_iCAPs,'subjectMaps',['iCAP_z_' num2str(param.K) '.nii']),'file');
            if isfield(param,'saveSubjectMaps') && param.saveSubjectMaps && ...
                    ~subjectsSaved
                load(fullfile(param.outDir_iCAPs,'IDX'));
                WriteInformation(fid,'Saving subject maps...');
                saveSubjectMaps(param,subject_labels,IDX,I_sig,final_mask);
            end
            
            %% saving tables with regions
            regTableExist=exist(fullfile(param.outDir_iCAPs,'iCAP_z_regions.txt'),'file');
            if isfield(param,'saveRegionTables') && param.saveRegionTables && ...
                    ~regTableExist
                WriteInformation(fid,'Saving region tables...');
                load(fullfile(param.outDir_iCAPs,'iCAPs_z'));
                saveRegionTables(param,iCAPs_z,final_mask);
            end
        end
    end
    
    % compute cluster stability based on consensus clustering (clustering
    % and consensus clustering have to be done already)
    if isfield(param,'computeClusterStability') && param.computeClusterStability
        WriteInformation(fid,'Getting cluster consensus...');
        % clustering for every K
        for iK=1:length(param.iCAPs_title_cell)
            param.iCAPs_title=param.iCAPs_title_cell{iK};
            param.K=param.K_vect(iK);
            WriteInformation(fid,['\nK=' num2str(param.K)]);
            param.outDir_iCAPs=fullfile(param.outDir_main,param.iCAPs_title);
            param.outDir_cons=fullfile(param.PathData,'iCAPs_results',...
                [param.data_title,'_',param.thresh_title],param.cons_title);
            % check files
            [~,Clustering_done,ConsensusClustering_done] = Check_iCAPs_Files(param.outDir_main,param.outDir_iCAPs,param.outDir_cons);
            if ~ConsensusClustering_done || ~Clustering_done
                WriteInformation(fid,'Run consensus clustering and clustering first, stability not computed!');
                continue
            end
            load(fullfile(param.outDir_iCAPs,'IDX'));
            load(fullfile(param.outDir_cons,['Consensus_' num2str(param.K)]));
            [iCAPs_consensus,iCAPs_nItems]=getClusterConsensus(IDX,Consensus);
            for iC=1:param.K
                WriteInformation(fid,['iCAP ' num2str(iC) ' (' num2str(iCAPs_nItems(iC)) ' frames) average consensus is ' num2str(iCAPs_consensus(iC))])
            end
            save(fullfile(param.outDir_iCAPs,'iCAPs_consensus'),'iCAPs_consensus');
            save(fullfile(param.outDir_iCAPs,'iCAPs_nItems'),'iCAPs_nItems');
        end
    end
end