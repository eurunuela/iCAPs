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
% - TC_OUT is the n_time_points x n_ret_voxels 2D matrix of outputs from
% the regularization step
% - param is the updated structure with relevant parameters for TA; added
% fields are 'LambdaTempFin' (vector with final regularization estimates
% for each voxel), 'NoiseEstimateFin' (final noise estimate for each voxel)
%
% Implemented by Isik Karahanoglu, 28.11.2011
function [TC_OUT,param] = TA_Temporal(TCN,param)
    
    % The output from the algorithm (time x voxels) is initialized as
    % a matrix of zeros
    TC_OUT = zeros(param.Dimension(4),param.NbrVoxels);
    
    % LambdaTemp contains the values of regularisation parameters for each
    % voxel; also set to zero for now
    param.LambdaTemp = zeros(param.NbrVoxels,1);
    
    % The noise estimation procedure uses a single scale
    % wavelet decomposition. Here Daubechies wavelets with 4 vanishing 
    % moments are used. The corresponding high pass filter is given by:

    g=[0    -0.12941    -0.22414     0.83652    -0.48296];
    g = g';


    % We iterate through all voxels to solve the problem
    for i=1:param.NbrVoxels,
        
        % Initialisation of the regularisation parameter for the
        % considered voxel, in two steps (a,b):

        % a. Wavelet decomposition of the time course of the voxel
        % of interest. 

        
        coef=cconv(TCN(:,i),g);


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

        [TC_OUT(:,i),paramOUT] = TA_Temporal_OneTimeCourse(TCN(:,i),i,param);

        
        % Takes the final estimates of regularisation parameter and
        % 'effective noise' for voxel i, and sroes them in the
        % param structure
        param.LambdaTempFin(i) = paramOUT.LambdasTempFin;
        param.NoiseEstimateFin(i) = paramOUT.NoiseEstimateFin;  
    end  
end

