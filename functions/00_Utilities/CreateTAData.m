%% This function creates a two-dimensional matrix containing the relevant
% data from total activation, from 4D input data and a mask determining the
% voxels to retain
%
% Inputs:
% - fData (X x Y x Z x T) is the functional data (4D)
% - mask is a 1D logical vector (prod(X,Y,Z)-long) with '1' if a voxel is
% part of the retained ones, and '0' else (i.e. out of the brain or NaN
% voxel)
% - param is a structure containing all the parameters relevant for total
% activation; here, we require the field 'Dimension' (4D vector with X, Y, 
% Z and T sizes), and we fill the fields 'IND' (1D indices of the elements
% to retain), 'VoxelIdx' (n_retained_voxels x 3 matrix with 3D coordinates 
% of the voxels to retain) and 'NbrVoxels' (scalar, how many voxels are
% kept)
% - fid points towards a log file where performed steps are registered
%
% Outputs:
% - fData_2D is a n_retained_voxels x n_time_points 2D matrix with the data
% used in TA
% - param is the updated parameters structure
function [fData_2D,param] = CreateTAData(fData,param,fid)

    fData_2D = nan(sum(param.mask),param.Dimension(4));

    for t = 1:param.Dimension(4)
        tmp = squeeze(fData(:,:,:,t));
        fData_2D(:,t) = tmp(param.mask);
    end
    
    % Finds the values of elements that are non-null both in the atlas and in
    % the data (1D vector)
    param.IND = find(param.mask==1);

    % Derives the 3D indices of those elements: param.VoxelIdx has dimensions n_el
    % x 3, with n_el the number of indices of non-null elements
    [param.VoxelIdx(:,1),param.VoxelIdx(:,2),param.VoxelIdx(:,3)] = ind2sub(param.Dimension(1:3),param.IND);

    WriteInformation(fid,['Keeping ',num2str(sum(param.mask)),' out of ',num2str(prod(param.Dimension(1:3))),' voxels for TA...']);
    
    param.NbrVoxels = size(fData_2D,1);
end