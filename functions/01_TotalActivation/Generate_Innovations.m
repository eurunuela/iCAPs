%% This function generates the innovation and activity inducing signals 
% from the activity related signals (output from TA)
%
% Inputs:
% - TC_OUT is a 2D matrix (n_time_points x n_ret_voxels) with the outputs
% from total activation
% - param contains the TA-relevant parameters; here, we need the fields
% 'Dimension' (size of X, Y, Z and T), 'NbrVoxels' (number of voxels
% retained for TA), 'f_Recons' (the filter that converts the
% activity-related signal into the activity-inducing one), 'f_Analyze' (the
% one with added derivation step)
%
% Outputs:
% - Innovation is a n_time_points x n_ret_voxels 2D matrix with the
% innovation signals
% - Activity_inducing is a n_time_points x n_ret_voxels 2D matrix with the
% activity-inducing signals
function [Innovation,Activity_inducing] = Generate_Innovations(TC_OUT,param)

    Activity_inducing = zeros(param.Dimension(4),param.NbrVoxels);
    Innovation = zeros(param.Dimension(4),param.NbrVoxels);

    % Each voxel time course is deconvolved using the reconstruction filter
    % (solely deconvolution)
    for i=1:param.NbrVoxels
        
        % Applies the reconstruction filter (only deconvolution)
        Activity_inducing(:,i) = filter_boundary(param.filter_reconstruct.num,param.filter_reconstruct.den,TC_OUT(:,i),'normal');
 
        % Applies the analysis filter (also encompasses the differentiation
        % step)
        Innovation(:,i) = filter_boundary(param.f_Analyze.num ,param.f_Analyze.den,TC_OUT(:,i),'normal');
        
        % Those are the ways originally present in the code: I have no idea
        % why they were used instead of the above option...
        %Innovation(:,i) = [0;diff((Activity_inducing(:,i)))];
        %Innovation(:,i) = cumsum([zeros(5,1); Innovation(6:end,i)]);
    end
    
end