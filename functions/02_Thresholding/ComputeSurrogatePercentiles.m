%% This function computes, for a set of innovation time courses, the
% percentile values at a given alpha-level (i.e. the threshold past which
% an innovation in real data is deemed significant)
%
% Inputs:
% - I (n_tp x n_ret_vox) is the matrix of surrogate data innovation signals
% - alpha is the level (in percentage) at which we want to threshold; it
% should be a 2-element vector (e.g. [1 99] for the 1st and 99th
% percentiles)
% - fid points towards the file where we want to write about the performed
% steps
%
% Outputs:
% - PC is a (2 x n_ret_vox) matrix with PC(1,:) the lower percentile and
% PC(2,:) the upper percentile values
function [PC] = ComputeSurrogatePercentiles(I,param,fid)

    % Computes percentiles for each time course (each column); PC is a 2 x
    % n_ret_vox matrix
    PC = prctile(I,param.alpha);
    
    % Writes in log-file
    WriteInformation(fid,['Computed voxel-wise percentiles: average values are ',...
        num2str(mean(PC(1,:))),' for bottom, and ',num2str(mean(PC(2,:))),...
        ' for top percentiles...']);
end