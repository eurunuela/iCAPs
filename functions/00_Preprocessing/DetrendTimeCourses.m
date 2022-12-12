%% This function detrends the time courses to consider for total activation
%
% Inputs:
% - TC is a n_ret_voxels x n_time_points 2D matrix with the data to detrend
% - param is the structure containing relevant TA parameters; here, we
% require the 'TR' (TR of the data), 'DCT_TS' (cutoff frequency for the DCT
% basis), 'Covariates' (list of covariates to include in the
% detrending/regression process; set to [] if none wished), and 'NbrVoxels'
% (number of retained voxels for TA)

% Younes, Jun 2019.
% Update: Allow for Normalization only.
% Younes, Dec 2019.
% Update: Save mean and Std maps from normalization.
function [TCN,STD_MAP] = DetrendTimeCourses(TC,param,fid)

    % TCN has n_vox x n_TP dimensions
    TCN = zeros(size(TC));
    STD_MAP = zeros(param.NbrVoxels,1);

    if  param.doDetrend

        for i=1:param.NbrVoxels;
            
            % Regresses out low frequency components + possible other
            % covariates
            TCN(i,:) = sol_dct(TC(i,:)',param.TR,param.DCT_TS,param.Covariates);
            
            % Normalization
            stdval = std(TCN(i,:));
            STD_MAP(i) = stdval;
            TCN(i,:) = TCN(i,:)./stdval;
        end

        WriteInformation(fid,['Detrending and normalizing the data with DCT = ',...
            num2str(param.DCT_TS),' [s] and ',num2str(size(param.Covariates,2)),...
            ' covariate(s)']);
    else 

        for i=1:param.NbrVoxels;      
            % Only Normalization
            stdval = std(TC(i,:));
            STD_MAP(i) = stdval;
            TCN(i,:) = TC(i,:)./stdval;
        end

        WriteInformation(fid,['Normalizing the data']);

    end

end