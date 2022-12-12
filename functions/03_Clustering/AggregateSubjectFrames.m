%% This function combines the frames from a set of subjects together into a
% full data matrix
%
% Inputs:
% - Paths is a cell array with each element a path towards the TA results 
% data from one subject
% - param contains the parameters of the TA routines
%
% - fid points towards the log file where applied operations are listed
%
% Outputs:
% - Data is a n_frames x n_final_ret_vox matrix with the frames to input to
% clustering
% - final_mask is the final mask across subjects considering only shared
% voxels
%
% 12.6.2018 - YF/DZ modified for different frame numbers in subjects

function [I_sig,final_mask,subject_labels,time_labels,AI,AI_subject_labels] = AggregateSubjectFrames(Paths,input_param,fid) 
    
    nSub=length(Paths);
    
    threshDir=fullfile(Paths,'Thresholding',input_param.thresh_title);
    taDir=fullfile(Paths,'TotalActivation');
    
%     clear Paths
    
    % Iterating through all subjects to consider...
    for iS = 1:nSub
        fprintf('%d ',iS)
        if ~exist(fullfile(threshDir{iS},'param.mat'),'file') ...
                || ~exist(fullfile(threshDir{iS},'SignInnov.mat'),'file') ...
                || ~exist(fullfile(taDir{iS},'Activity_inducing.mat'),'file')
            disp('......');
            error(['Cannot access data from subject ' Paths{iS} '; clustering process aborted...']);
        else
            % The mask is always the intersection of the current mask and the
            % new one (to remove elements that are not present everywhere)
            load(fullfile(threshDir{iS},'param.mat'),'param');
            % Data from the subject
            load(fullfile(threshDir{iS},'SignInnov.mat'));
            load(fullfile(taDir{iS},'Activity_inducing.mat'));
            
            if size(SignInnov,1)>size(SignInnov,2)
                warning(['SignInnov probably inverted - number of voxels: ' num2str(size(SignInnov,1)) ', number of frames ' num2str(size(SignInnov,2)) ', transposing data ...'])
                SignInnov=SignInnov';
            end
            if size(Activity_inducing,1)>size(Activity_inducing,2)
                warning(['Activity_inducing probably inverted - number of voxels: ' num2str(size(Activity_inducing,1)) ', number of frames ' num2str(size(Activity_inducing,2)) ', transposing data ...'])
                Activity_inducing=Activity_inducing';
            end
            
            % subject-specific data: nVox should be the same for all
            % subjects, nTP can be different and is only used for space
            % allocation here, nSignInnov is different for every subject
            nTP=length(param.mask_threshold2pos);
            nSignInnov=length(param.time_labels);
            nVox=length(param.mask);
            
            if iS == 1 % allocating memory in first iteration
                final_mask = ones(nVox,1);
                I_sig = sparse(2*nTP,nVox);
                AI = sparse(nTP,nVox);
                AI_subject_labels = zeros(nTP,1);
                subject_labels = zeros(2*nTP,1);
                time_labels = zeros(2*nTP,1);
                iI_sig=1;
                iAI=1;
            end
            
            % increase allocated size if necessary
            if size(I_sig,1)<iI_sig+nSignInnov-1
                I_sig=[I_sig;sparse(iI_sig+nSignInnov-1-size(I_sig,1),nVox)];
                subject_labels=[subject_labels;zeros(iI_sig+nSignInnov-1-size(I_sig,1),1)];
                time_labels=[time_labels;zeros(iI_sig+nSignInnov-1-size(I_sig,1),1)];
            end
            if size(AI,1)<iAI+nTP-1
                AI=[AI;sparse(iAI+nTP-1-size(AI,1),nVox)];
                AI_subject_labels=[AI_subject_labels;zeros(iAI+nTP-1-size(AI,1),1)];
            end
            
            % The final data matrix is updated
            % (n_frames_uptosubjectconsidered x n_vox)
            I_sig(iI_sig:iI_sig+nSignInnov-1,param.mask_nonan) = SignInnov;
            AI(iAI:iAI+nTP-1,param.mask) = Activity_inducing;
            
            % time and subject labels of significant innovations
            AI_subject_labels(iAI:iAI+nTP-1)=iS;
            subject_labels(iI_sig:iI_sig+nSignInnov-1) = iS;
            time_labels(iI_sig:iI_sig+nSignInnov-1) = param.time_labels;
            
            % The final mask is updated, removing elements if they are 0 in
            % the newly considered subject mask case
            final_mask = logical(final_mask & param.mask_nonan);
            
            % updating indices
            iI_sig=iI_sig+nSignInnov;
            iAI=iAI+nTP;
            
            clearvars param SignInnov Activity_inducing
        end
    end
    
    % remove remaining zeros from initialized matrices
    I_sig(iI_sig:end,:)=[];
    subject_labels(iI_sig:end,:)=[];
    time_labels(iI_sig:end,:)=[];
    AI(iAI:end,:)=[];
    
    % check if the common mask should be taken from an external file
    if isfield(input_param,'common_mask_file') && ~isempty(input_param.common_mask_file)
        final_mask=spm_read_vols(spm_vol(input_param.common_mask_file));
        final_mask=final_mask(:)>0.2;
    end
    
    % check if there is an additional mask file to apply on the
    % intersection of all GM masks
    if isfield(input_param,'extra_mask_file') && ~isempty(input_param.extra_mask_file)
        extra_mask_hdr=spm_vol(input_param.extra_mask_file);
        extra_mask=spm_read_vols(extra_mask_hdr);
        mapTo_hdr=spm_vol(fullfile(threshDir{iS},'SignInnov.nii,1'));
        extra_mask=mapVTV(extra_mask,extra_mask_hdr,mapTo_hdr);
        
        extra_mask=extra_mask(:)>0.2;
        final_mask=logical(final_mask & extra_mask);
    end
    
    % We can now mask the final data using the final mask
    I_sig = full(I_sig(:,final_mask));
    AI = full(AI(:,final_mask));
    
    param.final_mask = final_mask;
    
    % Writes down useful information in the log-file
    WriteInformation(fid,['\niCAPs clustering: Loaded the data from ',...
        num2str(nSub),' subjects for clustering...']);
    
    WriteInformation(fid,['There are in total ',num2str(size(I_sig,2)),...
        ' voxels kept for clustering...']);
    
    WriteInformation(fid,['There are in total ',num2str(size(I_sig,1)),...
        ' frames kept for clustering...']);
end