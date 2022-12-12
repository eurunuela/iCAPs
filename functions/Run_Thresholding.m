%% This function runs thresholding on innovations, according to the 
% parameters specified by the user
%
% Input: 
%   param - struct containing all necessary parameters to run TA
%       .PathData - path to data
%       .Subjects - cell array with list of subdirectories where each 
%           subject's fMRI data is stored; must contain one entry per 
%           subject to analyze.
%           This is where the TA folder will be created (or looked for) 
%           for each subject.
%       .n_subjects - number of subject to analyze
%       [.title] - possibility to define a title for the current
%           project (usefull if TA should be run for different 
%           parameters), default: current date
%       [.force_Thresholding] - if set to 1, Thresholding will be forced to 
%           run, even if already has been done
%       .alpha - Alpha-level at which to look for significance of 
%           innovation signal frames (first element is the percentile of 
%           the lower threshold - negative innovations - and second element 
%           the upper threshold one - positive innovations)
%       .f_voxels - Fraction of voxels from the ones entering total 
%           activation for a given subject that should show an innovation 
%           at the same time point, so that the corresponding frame is 
%           retained for iCAPs clustering
%       .thresh_title - Title used to create the folder where thresholding 
%           data will be saved
%       .threshold_minclussize - Number of neighbours that must also show 
%           an innovation for a voxel to be retained
%       .threshold_interconnectivity - Number of neighbours to consider in 
%           the process
%
% Output:
%       Cycles through the subjects and thresholds their data:
%       1) temporal threshold: takes the innovation signals generated from 
%       the surrogate data, builds a distribution from them, samples the 
%       X-th percentiles from the distribution at each voxel mask temporal 
%       frames from the real data with this
%       2) spatial threshold: mark frames where X percent of significant
%       voxels


function [] = Run_Thresholding(param)

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
    fid=fopen(fullfile(param.PathData,'TAlogs',['log_Thresholding_',param.title,'.txt']),'a+');
    WriteInformation(fid,['\nStarting the total activation/iCAPs tools for project entitled ',param.title]);
    
    % Checks that the path towards the data is correct
    if ~exist(param.PathData,'dir')
        WriteInformation(fid,'Incorrect path towards the data: execution stopped');
        error('The data folder that you specified does not exist ! Please check and restart running...');
    end
    
    
    
    WriteInformation(fid,'Entering the thresholding process...');
    
    
    % Loops over all the subjects to analyze
    for iS = 1:param.n_subjects
            
        % Path towards the data of the subject of interest
        SubjPath = fullfile(param.PathData,param.Subjects{iS});
        WriteInformation(fid,['Analyzing subject ',SubjPath,'...']);
        
        if ~exist(SubjPath,'dir')
            WriteInformation(fid,['Incorrect subject path ',SubjPath,': ignored subject']);
            continue
        end
        
        % Creates the folder that will contain total activation results
        % for the considered subject if it does not exist yet
        resultsPath=fullfile(SubjPath,'TA_results',param.title);
        
        
        % check whether thresholding was already executed for this subject
        [TA_real_done,TA_surrogate_done,Thresholding_done] = Check_TA_Files(resultsPath,param.thresh_title);
        if ~TA_real_done || TA_surrogate_done~=1
           WriteInformation(fid,['No total activation results, run TA routine first! skipping...\n']);
           continue;
        end
        
        if isfield(param,'force_Thresholding') && param.force_Thresholding
            Thresholding_done=0;
        end
        
        if ~Thresholding_done
        
            % loading data
            WriteInformation(fid,['Loading total activation data...']);
            load(fullfile(resultsPath,'TotalActivation','Innovation'));
            load(fullfile(resultsPath,'Surrogate','Innovation_surrogate'));
            ptmp=load(fullfile(resultsPath,'TotalActivation','param'));

            param.fHeader=ptmp.param.fHeader;
            param.mask=ptmp.param.mask;
            param.Dimension=ptmp.param.Dimension;
            if isfield(ptmp.param,'TemporalMask')
                param.TemporalMask=ptmp.param.TemporalMask;
            end
            clear ptmp

            WriteInformation(fid,['Removing NaNs from innovations...']);
            [Innovation_surrogate,param.mask_nonan,param.mask2_nonan] = RemoveNan(Innovation_surrogate,param,fid);
            
            WriteInformation(fid,['Computing percentiles...']);
            param.PC = ComputeSurrogatePercentiles(Innovation_surrogate,param,fid);

            WriteInformation(fid,['Selecting significant innovation frames...']);
            [SignInnov,param] = SelectSignificantFrames(Innovation(:,param.mask2_nonan),param,fid);

            WriteInformation(fid,['Saving...']);
            save4Dnii(fullfile(resultsPath,'Thresholding'),param.thresh_title,'SignInnov',SignInnov',param.fHeader.fname,param.mask_nonan,param.Dimension);
            save4Dnii(fullfile(resultsPath,'Thresholding'),param.thresh_title,'mask_nonan',param.mask_nonan,param.fHeader.fname,param.mask,param.Dimension);
            save(fullfile(resultsPath,'Thresholding',param.thresh_title,'SignInnov'),'SignInnov','-v7.3');
            save(fullfile(resultsPath,'Thresholding',param.thresh_title,'param'),'param','-v7.3');
            
            WriteInformation(fid,['Finished running thresholding for subject ',SubjPath,'...\n']);
            
        elseif Thresholding_done==1
             WriteInformation(fid,'Thresholding already done, skipping...\n');
        end
            
        
        % Resets the parameters to what they were at the start of the loop
        % (before any subject-specific change could have been made
        clear param
        param = param_CI;
    end
end