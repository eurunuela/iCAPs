%% This function performs the double thresholding process: (1) assess the
% moments of an innovation time course with significant excursions, and (2)
% assess the time points with sufficient such excursions brain-wise
%
% Inputs:
% - Innovation (n_tp x n_ret_vox) contains the innovation time courses for
% the considered subject (real data)
% - PC is a 2 x n_ret_vox matrix with percentile values for each voxel, as
% computed on surrogate data
% - f_voxels is the fraction of voxels that must be showing a significant
% excursion at the same time to retain a time point for the iCAPs
% clustering
% - param is the structure containing relevant TA parameters; it will be
% updated with additional entries
%
% Outputs:
% - SignInnov is a n_significant_excursions x n_ret_vox matrix containing
% the frames which showed significant innovation
function [SignInnov,param] = SelectSignificantFrames(Innovation,param,fid)

    % Will be a matrix specifying the time points with significant
    % innovations for each voxel (-1 for negative excursion, 0 for no
    % excursion, +1 for positive excursion)
    param.mask_threshold1 = ThresholdTimeCourse(Innovation,param.PC);
    
%     WriteInformation(fid,['Temporal thresholding: there are on average ',...
%         num2str(mean(sum(abs(param.mask_threshold1),1))),...
%         ' out of ',num2str(size(Innovation,1)),' time points kept...']);
%     
%     WriteInformation(fid,['Positive time points: ',...
%         num2str(mean(sum(param.mask_threshold1==1,1)))]);
%     
%     WriteInformation(fid,['Negative time points: ',...
%         num2str(mean(sum(param.mask_threshold1==-1,1)))]);
        
    % Computes how many voxels must show significance at once
    n_voxels = floor(size(Innovation,2)*param.f_voxels);
    
    WriteInformation(fid,['Spatial thresholding: there are ',...
        num2str(n_voxels),...
        ' voxels that must show significance together...']);
    
    % Mask only retaining voxels/time points with POSITIVE innovations (+1
    % elements; all other elements set to 0)
    mask_threshold1pos=param.mask_threshold1;
    mask_threshold1pos(mask_threshold1pos<0)=0;

    % Same for negative elements; the mask is filled with only -1 (negative
    % innovations) or 0 (nothing or positive innovation)
    mask_threshold1neg=param.mask_threshold1;
    mask_threshold1neg(mask_threshold1neg>0)=0;
        
    % Further processes the masks so that the too small voxel islands of
    % innovation (less than 6 neighbours in 3D) are removed
    mask_threshold1bpos = check_interconnectedness(mask_threshold1pos,param);
    mask_threshold1bneg = check_interconnectedness(mask_threshold1neg,param);

    % Generates a n_timepoints x 1 vector indexing what frames to retrieve
    % because it is a positive (pos) or negative (neg) innovation
    param.mask_threshold2pos = ThresholdWholeBrain(mask_threshold1bpos,n_voxels);
    param.mask_threshold2neg = ThresholdWholeBrain(mask_threshold1bneg,n_voxels);

    WriteInformation(fid,['There are ',...
        num2str(sum(param.mask_threshold2pos)),...
        ' out of ',num2str(size(Innovation,1)),' frames selected for positive innovation...']);
    
    WriteInformation(fid,['There are ',...
        num2str(sum(param.mask_threshold2neg)),...
        ' out of ',num2str(size(Innovation,1)),' frames selected for negative innovation...']);
    
    % Stores the frames showing positive innovation (which have had their
    % negative elements removed)
    SignInnovPos = Innovation(param.mask_threshold2pos,:);
    SignInnovPos(SignInnovPos<0) = 0;
    
    % Stores the frames showing negative innovation (which have had their
    % sign flipped, and their positive inovation elements removed)
    SignInnovNeg = -Innovation(param.mask_threshold2neg,:);
    SignInnovNeg(SignInnovNeg<0) = 0;
    
    % Final data matrix
    SignInnov=[SignInnovPos;SignInnovNeg];
    
    param.time_labels=[find(param.mask_threshold2pos==1)',find(param.mask_threshold2neg==1)'];
    
end