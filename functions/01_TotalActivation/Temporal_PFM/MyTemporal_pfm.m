%% This function performs the temporal regularization part of total 
% activation, voxel after voxel
%
% Inputs:
% - TCN is the n_time_points x n_ret_voxels 2D matrix of data input to the
% regularization
% - param is a structure containing all TA-relevant parameters; here, we
% will need the fields 'Dimension' (X, Y, Z and T sizes), 'NbrVoxels'
% (number of voxels to consider for regularization), 'LambdaTempCoef' (used
% to compute regularization coefficients)
%
% Outputs:
% - Activity_related is the n_time_points x n_ret_voxels 2D matrix of outputs from
% the regularization step
% - param is the updated structure with relevant parameters for TA; added
% fields are 'LambdaTempFin' (vector with final regularization estimates
% for each voxel), 'NoiseEstimateFin' (final noise estimate for each voxel)
%
% Implemented by Eneko Uruñuela, 13.12.2022
function [Activity_related,Activity_inducing,innovation, paramOUT] = MyTemporal_pfm(TCN, param)

    % The output from the algorithm (time x voxels) is initialized as
    % a matrix of zeros
    Activity_related = zeros(param.Dimension(4),param.NbrVoxels);
    Activity_inducing = zeros(param.Dimension(4),param.NbrVoxels);
    innovation = zeros(param.Dimension(4),param.NbrVoxels);
    % LambdaTemp contains the values of regularisation parameters for each
    % voxel; also set to zero for now
    param.LambdaTemp = zeros(param.NbrVoxels,1);

    % Generate HRF
    param.HRF = GenerateHRF(param.TR, param.Dimension(4), param.block, param.custom);

    % The necessary HRF matrices are computed
    param.X_tilde = param.HRF;
    param.X_tilde_trans = param.X_tilde';
    param.X_tilde_tt = param.X_tilde_trans*param.X_tilde;
    param.c_ist = 1/abs(eigs(param.X_tilde_tt,1));

    % We loop over all voxels
    for i = 1:param.NbrVoxels

        % If LambdaPFM is "mad", then we use the median absolute deviation
        % to estimate the noise level
        
        if strcmp(param.LambdaPFM,'mad')

            % a. Wavelet decomposition of the time course of the voxel
            % of interest
            [coef,len] = wavedec(TCN(:,i),1,'db3');

            % b. The output from the wavelet decomposition is a separation
            % of the signal into an 'approximate' and a 'detailed'
            % part, with len(1) the number of 'approximate' coefficients,
            % Those are removed from the 'coef' vector to keep only
            % noise-like coefficients
            coef(1:len(1)) = [];

            % c. Median absolute deviation (sum of absolute valued
            % distances of coefficients from the mean)
            % From the Matlab page on wavelet denoising:
            % "The median absolute deviation of the coefficients is a
            % robust estimate of noise."
            % This is thus going to be our estimate of the noise level
            % for the considered voxel time course
            param.LambdaTemp(i) = mad(coef,1)*param.LambdaTempCoef;

             % Now that we have estimated our initial lambda 
            % (regularisation parameter), Temporal_TA performs the 
            % computations themselves for the considered time course
            % TCN(i,:) of voxel i

            [Activity_related(:,i),Activity_inducing(:,i),innovation(:,i),paramOUT] = PFM_Temporal_OneTimeCourse(TCN(:,i),i,param);
        
        % If LambdaPFM is "aic" or "bic", use PFM with LARS
        else strcmp(param.LambdaPFM,'aic') || strcmp(param.LambdaPFM,'bic')
            % Error message saying that this is not implemented yet
            error('PFM with LARS is not implemented yet')
        end


    end
end