%% This function runs regression to obtain iCAPs time courses
%  regression can either be done 'unconstrained' or with 'transient-based'
%  constraints
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
%       [.force_Regression] - if set to 1, regression will be forced to 
%           run, even if already has been done
%       [.thresh_title] - information to read thresholding data and save
%           regression results, if not specified, the fields 'param.alpha'
%           and 'param.f_voxels' need to exist
%       [.data_title] - information for saving of regression results
%       [.iCAPs_title] - string or cellstring with all the subfolders to
%           create for saving regression results
% 
%     * Regression information:
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
% Output:


function [] = Run_Regression(param)

    % Date and time when the routines are called
    param.date = strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
    if ~isfield(param,'title') || isempty(param.title)
        param.title=param.date;
    end
    
    % Creates and opens a log-file that will contain all information related to
    % what is done to the data
    if ~exist(fullfile(param.PathData,'TAlogs'),'dir'); mkdir(fullfile(param.PathData,'TAlogs'));end;
    fid=fopen(fullfile(param.PathData,'TAlogs',['log_Regression_',param.title,'.txt']),'a+');
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
    
    % main folder according to data+thresholding, aggregated frames are saved here
    param.outDir_main=fullfile(param.PathData,'iCAPs_results',[param.data_title,'_',param.thresh_title]);
    
    % iCAPs title (subfolder in which clustering results are saved)
    if ~isfield(param,'iCAPs_title') || isempty(param.iCAPs_title)
        for nK=1:length(param.K)
            param.iCAPs_title{nK} = ['K_',num2str(param.K(nK)),'_Dist_',...
                    param.DistType,'_Folds_',num2str(param.n_folds)];
        end
        if length(param.K)==1
            param.iCAPs_title=param.iCAPs_title{1};
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
    
    % main procedure for regression:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(param,'doRegression') || param.doRegression
        WriteInformation(fid,'Entering the time course retrieval process...');
        % clustering for every K
        for iK=1:length(param.iCAPs_title_cell)

            param.iCAPs_title=param.iCAPs_title_cell{iK};
            param.K=param.K_vect(iK);

            WriteInformation(fid,['K = ' num2str(param.K) ', iCAPs title: ' param.iCAPs_title]);
            
            % clustering subfolder to save clustering-specific information
            param.outDir_iCAPs=fullfile(param.outDir_main,param.iCAPs_title);
            
            % check if clustering has been done
            [Aggregating_done,Clustering_done,~,~] = Check_iCAPs_Files(param.outDir_main,param.outDir_iCAPs);

            if ~Aggregating_done
                error('Significant innvation frames have not been aggregated yet, run clustering routine first!')
            end
            if ~Clustering_done
                error('Clustering has not been done yet, run clustering routine first!')
            end
            
                
            switch param.regType
                case 'unconstrained'
                    param.outDir_reg=fullfile(param.outDir_iCAPs,'TCs_unconstrained');
                    if ~exist(param.outDir_reg,'dir');mkdir(param.outDir_reg);end
                    [~,~,~,~,Regression_unc_done] = Check_iCAPs_Files([],[],[],param.outDir_reg);

                    if isfield(param,'force_Regression') && param.force_Regression
                        Regression_unc_done=0;
                    end
                    if ~Regression_unc_done
                        WriteInformation(fid,'Loading aggregated data and iCAPs...');
                        load(fullfile(param.outDir_main,'AI.mat'));
                        load(fullfile(param.outDir_main,'AI_subject_labels.mat'));
                        load(fullfile(param.outDir_iCAPs,'iCAPs'));
                        
                        WriteInformation(fid,'Generating iCAPs Time Courses...');
                        [TC_unc,TC_unc_stats] = GenerateTimeCourses(AI',AI_subject_labels,iCAPs,param);
                        
                        save(fullfile(param.outDir_reg,'TC_unc.mat'),'TC_unc');
                        save(fullfile(param.outDir_reg,'TC_unc_stats.mat'),'TC_unc_stats');
                        save(fullfile(param.outDir_reg,'param.mat'),'param');
                        
                        load(fullfile(param.outDir_main,'AI_subject_labels'));
                        load(fullfile(param.outDir_main,'subject_labels'));
                        load(fullfile(param.outDir_iCAPs,'IDX'));
                        clusteringResults.AI_subject_labels=AI_subject_labels;
                        clusteringResults.subject_labels=subject_labels;
                        clusteringResults.IDX=IDX;
                        
                        tempChar_unc =computeTemporalCharacteristics(TC_unc,clusteringResults,param);
                        save(fullfile(param.outDir_reg,'tempChar_unc.mat'),'tempChar_unc');
                        
                    else
                        WriteInformation(fid,'Unconstrained regression has already been done, skipping...');
                    end
                case 'transient-informed'
                    param.outDir_reg=fullfile(param.outDir_iCAPs,...
                        strrep(['TCs_' num2str(param.softClusterThres(1)) '_' ...
                        num2str(mean(diff(param.softClusterThres))) '_' ...
                        num2str(param.softClusterThres(end))],'.','DOT'));
                    if ~exist(param.outDir_reg,'dir');mkdir(param.outDir_reg);end
                    [~,~,~,Regression_ti_done,~] = Check_iCAPs_Files([],[],[],param.outDir_reg);
                    
                    if isfield(param,'force_Regression') && param.force_Regression
                        Regression_ti_done=0;
                    end
                    if ~Regression_ti_done
                        WriteInformation(fid,'Loading aggregated data and iCAPs...');
                        clusteringResults=Load_ClusteringResults(param);% the first input variable is not used anymore, since I compute the path from param

                        for iT=1:length(param.softClusterThres)
                            WriteInformation(fid,['Soft assignment factor: ' num2str(param.softClusterThres(iT))]);
                            param_tmp=param;
                            param_tmp.softClusterThres=param.softClusterThres(iT);
                            [TC{iT},TC_stats{iT}]=GenerateTimeCoursesWeighted(clusteringResults,param_tmp);
                            tempChar{iT} =computeTemporalCharacteristics(TC{iT},clusteringResults,param_tmp);
                        end
                        save(fullfile(param.outDir_reg,'TC.mat'),'TC');
                        save(fullfile(param.outDir_reg,'TC_stats.mat'),'TC_stats');
                        save(fullfile(param.outDir_reg,'param.mat'),'param');
                        save(fullfile(param.outDir_reg,'tempChar.mat'),'tempChar');
                        
                        if isfield(param,'evalAmplitudeCorrs') && param.evalAmplitudeCorrs
                            evaluateSoftClusterThres_corrs(clusteringResults,TC,param,fid)
                        end
                        
                        % add a function to plot BIC and AIC (and maybe the
                        % correlation with measured amplitudes?)
                        best_id_bic=evaluateSoftClusterThres(TC_stats,param,fid);
                        TC=TC{best_id_bic};
                        TC_stats=TC_stats{best_id_bic};
                        tempChar=tempChar{best_id_bic};
                        best_val=strrep(num2str(param.softClusterThres(best_id_bic)),'.','DOT');
                        save(fullfile(param.outDir_iCAPs,['TC_' best_val '.mat']),'TC');
                        save(fullfile(param.outDir_iCAPs,['TC_stats_' best_val '.mat']),'TC_stats');
                        save(fullfile(param.outDir_iCAPs,['tempChar_' best_val '.mat']),'tempChar');
                        
                    else
                        WriteInformation(fid,'Transient-informed regression has already been done, skipping...');
                    end
                        
            end
                
                
        end
    end
end