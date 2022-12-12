%% This function takes outputs from total activation, and adjusts the mask
% to remove NaN time courses. Those are due to voxels from the structural
% data that were not in the functional acquisition volume; a flat trace at
% zero leads to NaN outcomes out of the TA pipeline
%
% Inputs:
% - Data is a n_TP x n_ret_vox matrix of voxels
% - mask is the n_vox x 1 vector specifying what voxels have been inputs to
% the TA analysis
%
% Outputs:
% - Data2 is a n_TP x n_ret_vox_2 matrix with the NaN traces removed
% - mask2 is a new mask with the NaN time courses removed
function [Innovation2,mask_out,mask2] = RemoveNan(Innovation,param,fid)
    
    if size(Innovation,2)<size(Innovation,1)
        warning(['data probably inverted - number of voxels: ' num2str(size(Innovation,2)) ', number of frames ' num2str(size(Innovation,1)) ', transposing data ...'])
        Innovation=Innovation';
    end
    
    % If the time course has no NaN values, we want to include it in
    % the data and the mask is filled appropriately
    mask2=sum(isnan(Innovation))==0;
    Innovation2=Innovation(:,mask2);
    
    % Fills the output mask that contains information for all voxels
    mask_out=param.mask;
    mask_out(param.mask)=mask2;
    
    WriteInformation(fid,['There are ',num2str(sum(mask2==0)),...
        ' voxels that have been removed because of being NaN after TA...']);
end