function [TC_OUT, param] = PFM_temporal(TCN, param)

    % The output from the algorithm (time x voxels) is initialized as
    % a matrix of zeros
    TC_OUT = zeros(param.Dimension(4),param.NbrVoxels);

    % LambdaTemp contains the values of regularisation parameters for each
    % voxel; also set to zero for now
    param.LambdaTemp = zeros(param.NbrVoxels,1);

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
        
        % If LambdaPFM is "aic" or "bic", use PFM with LARS
        elseif strcmp(param.LambdaPFM,'aic') || strcmp(param.LambdaPFM,'bic')
            
            % We use the PFM_LARS function to estimate the
            % time course of the voxel of interest
            [TC_OUT(:,i), param.LambdaTemp(i)] = PFM_LARS(TCN(:,i),param);

        % Else, use the universal threshold by default
        else
            
            % We use the PFM_universal function to estimate the
            % time course of the voxel of interest
            [TC_OUT(:,i), param.LambdaTemp(i)] = PFM_universal(TCN(:,i),param);
            
        end


    end

end