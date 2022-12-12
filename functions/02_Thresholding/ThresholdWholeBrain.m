%% This function determines the time points that exhibit whole-brain
% significant excursions in terms of innovation signal
%
% Inputs:
% - mask1 is a n_TP x n_ret_vox matrix containing -1 or 1 for significant
% negative/positive excursions at a given time for a given voxel, and 0 for
% no significant excursion
% - n_voxels is the number of voxels that must show significant excursions
% at the same time to conclude significance
%
% Outputs:
% - mask2 is a n_TP-long vector depicting the moments with significant
% excursions
function [mask2] = ThresholdWholeBrain(mask1,n_voxels)
    
    mask2 = zeros(size(mask1,1),1);
    mask2(sum(abs(mask1),2) >= n_voxels) = 1;
    mask2 = logical(mask2);
end