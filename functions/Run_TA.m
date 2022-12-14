%% This function runs total activation, according to the 
% parameters specified by the user
%
% Input: 
%   param - struct containing all necessary parameters to run TA
%       .PathData - path to data
%       .TR - TR of the fMRI data
%       .Subjects - cell array with list of subdirectories where each 
%           subject's fMRI data is stored; must contain one entry per 
%           subject to analyze.
%           This is where the TA folder will be created (or looked for) 
%           for each subject.
%       .n_subjects - number of subject to analyze
%       [.title] - possibility to define a title for the current
%           project (usefull if TA should be run for different 
%           parameters), default: current date
%       [.force_TA_on_real] - if set to 1, TA will be forced to run,
%           even if already has been done
%       [.force_TA_on_surrogate] - if set to 1, TA on surrogate data
%           will be forced to run, even if already has been done
%       .Folder_functional - name of the functional folder; [] if directly 
%           lying in Path
%       .TA_func_prefix - string with the prefix for functional data to read
%       .Folder_GM - name of the folder with the probabilistic gray matter map
%       .TA_GM_prefix - string with the prefix of the probabilistic map to read
%       .T_gm - threshold probability for creating GM mask
%       [.is_morpho] - if 1, morphological operations (opening and closure)
%           will be run on the GM mask to remove wholes
%       [.n_morpho_voxels] - required if 'is_morpho' is set, size of
%           morphological opterators
%       .skipped_scans - number of fMRI scans to skip at the beginning
%       .doDetrend - select if detrending should be done
%       .doNormalize - select if normalizing without detrensing
%       [.DCT_TS] - required if 'doDetrend' is set, cut-off period for the 
%           DCT basis
%       [.Covariates] - required if 'doDetrend' is set, covariates to add
%       .doScrubbing - select if motion censoring should be done
%       [.Folder_motion] -  required if 'doScrubbing' is set, folder with 
%           motion data from SPM realignment
%       [.TA_mot_prefix] - required if 'doScrubbing' is set, prefix of
%           motion data text file
%       [.skipped_scans_motionfile] - Number of lines to ignore at the 
%           beginning of the motion file,if empty or not set, this will be 
%           equal to param.skipped_scans
%       [.FD_method] - required if 'doScrubbing' is set, for now, only
%           'Power' is implemented
%       [.FD_threshold] - required if 'doScrubbing is set, scrubbing
%           threshold in mm
%       [.interType] - interpolation method (see also interp1), default is
%           'spline'
%
% Output:
%       creates a folder TA_results/<title> in each subject's folder and
%       saves results from total activation routine in subfolders:
%       - inputData: data after preprocessing
%       - TotalActivation: results after running total activation on real
%       - Surrogate: results after running total activation on surrogate data


function [] = Run_TA(param)

    % Date and time when the routines are called
    param.date = strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
    if ~isfield(param,'title') || isempty(param.title)
        param.title=param.date;
    end
    
    % Will contain the parameters initially entered by the user, prior to any
    % change within the loop
    param_CI = param;
    
    % Creates and opens a log-file that will contain all information related to
    % what is done to the data
    if ~exist(fullfile(param.PathData,'TAlogs'),'dir'); mkdir(fullfile(param.PathData,'TAlogs'));end;
    fid=fopen(fullfile(param.PathData,'TAlogs',['log_TA_',param.title,'.txt']),'a+');
        if fid==-1
      error('Cannot open file for writing ! Please check permission: %s', file);
        end
    WriteInformation(fid,['Starting the total activation/iCAPs tools for project entitled ',param.title]);
    
    % Checks that the path towards the data is correct
    if ~exist(param.PathData,'dir')
        WriteInformation(fid,'Incorrect path towards the data: execution stopped');
        error('The data folder that you specified does not exist ! Please check and restart running...');
    end

    
    % If the user wants to deploy the total activation parts of the script,
    % we enter a loop over the chosen subjects...
    WriteInformation(fid,'Entering the total activation part of the routines...');
        
    for i_TA = 1:param.n_subjects

        % Gathering subject information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Path towards the data of the subject of interest
        SubjPath_TA = fullfile(param.PathData,param.Subjects{i_TA});
        WriteInformation(fid,['Analyzing subject ',SubjPath_TA,'...']);
        
        if ~exist(SubjPath_TA,'dir')
            WriteInformation(fid,['Incorrect subject path ',SubjPath_TA,': ignored subject']);
            continue
        end
        
        % Creates the folder that will contain total activation results
        % for the considered subject if it does not exist yet
        resultsPath=fullfile(SubjPath_TA,'TA_results',param.title);
        if ~exist(resultsPath)
           mkdir(resultsPath);  
        end
        
        % check whether TA was already executed for this subject
        [TA_real_done,TA_surrogate_done] = Check_TA_Files(resultsPath);
        
        if isfield(param,'force_TA_on_real') && param.force_TA_on_real
            TA_real_done=0;
        end
        if isfield(param,'force_TA_on_surrogate') && param.force_TA_on_surrogate
            TA_surrogate_done=0;
        end

        if ~TA_real_done || TA_surrogate_done~=1
            % Reading subject data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Reads the data required for TA, i.e. the functional volumes and the
            % probabilistic GM map
            [fData,pData,fHeader,pHeader] = ReadTAData(SubjPath_TA,i_TA,param,fid);
            param.GM_map=pData;
            
            clear pData pHeader
            
            % Preprocessing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Dicards the first volumes from the functional data
            fData(:,:,:,1:param.skipped_scans)=[];
            WriteInformation(fid,['Discarding the first ',num2str(param.skipped_scans),' volumes from the time courses...']);

            % Dimension of the data (X x Y x Z x T)
            param.Dimension = size(fData);
            param.fHeader = fHeader;
            
%             % saving input fMRI data as nifti file
%             WriteInformation(fid,'Saving original fMRI 4D input (fData)...');
%             save4Dnii(resultsPath,'inputData','fData',fData,param.fHeader);
             
            % Creates the mask that will be used for the analysis: we want to keep
            % only the brain information
            [param.mask,param.mask_3D] = CreateTAMask(param,fid);
            
            % Creates the 2D data used for most of total activation
            [TC,param] = CreateTAData(fData,param,fid);
            
            % saving mask as nifti file
            WriteInformation(fid,'Saving mask 4D input...');
            save4Dnii(resultsPath,'inputData','mask',param.mask_3D,param.fHeader.fname);
            save4Dnii(resultsPath,'inputData','fData',fData,param.fHeader.fname,param.mask,param.Dimension);
            
            clear fData fHeader
            
            % Motion Analysis
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if param.doScrubbing
                param.TemporalMask = AssessMotion(SubjPath_TA,i_TA,param,fid);

                % This if condition is entered if there is at least one value in
                % TemporalMask that is 0 (i.e. if there is at least one frame for which we
                % must interpolate)
                if ~all(param.TemporalMask)
                    [TC, param.TemporalMask] = InterpolateTimeCourses(TC,param.TemporalMask,param,fid);
                else
                    WriteInformation(fid,['No interpolation done for ',SubjPath_TA,'...']);
                end
            end
            
            % Detrending
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if param.doDetrend || param.doNormalize
                % Detrended time courses (n_vox x n_TP)
                [TC,STD_MAP] = DetrendTimeCourses(TC,param,fid);
            end
            
            % Update time-course length after interpolation
            % because of excluded last frames
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            param.Dimension(4)= size(TC,2);
            
            
            % saving preprocessed input data as 4D nifti file
            WriteInformation(fid,'Saving preprocessed fMRI 4D input (TC)...');
            save4Dnii(resultsPath,'inputData','TC',TC,param.fHeader.fname,param.mask,param.Dimension);
            save4Dnii(resultsPath,'inputData','STD_MAP',STD_MAP,param.fHeader.fname,param.mask,param.Dimension);
            
            
            
            % Total activation itself
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The 'analyze' and 'reconstruct' filter cases are
            % different: 'analyze' has an added zero, which means, an additional
            % derivative (probably reflects the fact that we are working with
            % sparsity imposed at the level of the innovation signal, the
            % derivative of the piece-wise constant neural activity)
            
            % Creates the filters required: the one that 'deconvolves the BOLD
            % signal and derivates it' (analyze), and the one that 'deconvolves
            % the the BOLD-like signal into neural activity' (reconstruct). The
            % matlab variables contain the non-null values of the filter 
            % coefficients, from sample f[n] = f[0]
            param = hrf_filters(param);

            % Make the graph from gray matter voxels for spatial regularization
         
            param = make_graph(param);
            param.nu = 1;                      %gradient step for spatial reg


            % The param vector is updated within the total activation scheme; I
            % want to give exactly the same input for the surrogate and the real
            % data cases, so I save the state of the param vector prior to TA
            param_tmp = param;
            
            % Runs total activation first for the real data, and second for the
            % surrogate data
            
            % Real data
            if TA_real_done==0
                tic;
                % running total activation
                [Activity_related,param] = RunTotalActivation(TC',param);
                param.time_real = toc;
                WriteInformation(fid,['It took ',num2str(param.time_real),' seconds to run total activation on real data...']);
                
                % TA has been run, so now we can derive the activity-inducing and 
                % innovation signals from the activity related signal
                if param.use_PFM == 1
                    Innovation = param.Innovation;
                    Activity_inducing = param.Activity_inducing;
                else
                    [Innovation,Activity_inducing] = Generate_Innovations(Activity_related,param);
                end
                
                % Saving results from real data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                WriteInformation(fid,'Saving total activation results (mat and nifti)...');
                save4Dnii(resultsPath,'TotalActivation','Activity_inducing',Activity_inducing',param.fHeader.fname,param.mask,param.Dimension);
                save4Dnii(resultsPath,'TotalActivation','Activity_related',Activity_related',param.fHeader.fname,param.mask,param.Dimension);
                save4Dnii(resultsPath,'TotalActivation','Innovation',Innovation',param.fHeader.fname,param.mask,param.Dimension);
                save(fullfile(resultsPath,'TotalActivation','Activity_inducing'),'Activity_inducing','-v7.3');
                save(fullfile(resultsPath,'TotalActivation','Activity_related'),'Activity_related','-v7.3');
                save(fullfile(resultsPath,'TotalActivation','Innovation'),'Innovation','-v7.3');
                save(fullfile(resultsPath,'TotalActivation','param'),'param','-v7.3');

                % clear memory
                clear Activity_related Activity_inducing Innovation STD_MAP
                    
            elseif TA_real_done==1
                WriteInformation(fid,'Total activation on real data already computed, skipping...');
            end
            
            % surrogate data
            if TA_surrogate_done~=1
                % Surrogate data generation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Surrogate = GenerateSurrogate(TC,SubjPath_TA,param,fid);
                clear TC
                Surrogate=Surrogate';
                
                % saving surrogate data
                save4Dnii(resultsPath,'Surrogate','Surrogate',Surrogate',param.fHeader.fname,param.mask,param.Dimension);
                save(fullfile(resultsPath,'Surrogate','Surrogate'),'Surrogate','-v7.3');
                
                % Running TA for Surrogates
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tic;
                Activity_related_surrogate = RunTotalActivation(Surrogate,param_tmp);
                param.time_surrogate = toc;
                WriteInformation(fid,['It took ',num2str(param.time_surrogate),' seconds to run total activation on surrogate data...']);
                
                % TA has been run, so now we can derive the activity-inducing and 
                % innovation signals from the activity related signal
                if param.use_PFM == 1
                    Innovation_surrogate = param.Innovation;
                    Activity_inducing_surrogate = param.Activity_inducing;
                else
                    [Innovation_surrogate,Activity_inducing_surrogate] = Generate_Innovations(Activity_related_surrogate,param);
                end

                % Saving results from Surrogates
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                WriteInformation(fid,'Saving total activation results of surrogate data (mat and nifti)...');
                save4Dnii(resultsPath,'Surrogate','Activity_inducing_surrogate',Activity_inducing_surrogate',param.fHeader.fname,param.mask,param.Dimension);
                save4Dnii(resultsPath,'Surrogate','Activity_related_surrogate',Activity_related_surrogate',param.fHeader.fname,param.mask,param.Dimension);
                save4Dnii(resultsPath,'Surrogate','Innovation_surrogate',Innovation_surrogate',param.fHeader.fname,param.mask,param.Dimension);
                save(fullfile(resultsPath,'Surrogate','Activity_inducing_surrogate'),'Activity_inducing_surrogate','-v7.3');
                save(fullfile(resultsPath,'Surrogate','Activity_related_surrogate'),'Activity_related_surrogate','-v7.3');
                save(fullfile(resultsPath,'Surrogate','Innovation_surrogate'),'Innovation_surrogate','-v7.3');
                save(fullfile(resultsPath,'Surrogate','param'),'param','-v7.3');

                % Clear memory
                clear Activity_related_surrogate Activity_inducing_surrogate Innovation_surrogate Surrogate
                
            elseif TA_surrogate_done==1
                 WriteInformation(fid,'Total activation on surrogate data already computed, skipping...');
            end
            
            WriteInformation(fid,['Finished running total activation for subject ',SubjPath_TA,'...']);
        end
        
        % Resets the parameters to what they were at the start of the loop
        % (before any subject-specific change could have been made
        clear param
        param = param_CI;
        
    end
end