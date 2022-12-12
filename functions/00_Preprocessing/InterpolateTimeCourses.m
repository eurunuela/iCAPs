%% This function interpolates the fMRI volumes that are scrubbed out from
% the data, using spline interpolation
%
% Inputs:
% - fData is the functional data to consider (2D matrix, n_ret_vox x n_tp)
% - TemporalMask is a 1D logical vector with '0' for the time points where
% scrubbing is desired
%
% Outputs:
% - TC is the output 2D matrix (n_ret_vox x n_tp) following interpolation
function [TC,TemporalMask] = InterpolateTimeCourses(fData,TemporalMask,param,fid)
    
    if ~isfield(param,'interType') || isempty(param.interType)
        param.interType='spline';
    end
    nTP=length(TemporalMask);
    % if the last (or first) frame is excluded, don't include the last (or first) frames in the
    % interpolation
    indd=find(~TemporalMask);
    if any(indd==1)
        diffMask=diff(indd);
        excludeTP=1:indd(find(diffMask>2,1,'first')); % I choose 2 as a threshold here, in case one scan remains between excluded scans, here I will skip this one as well because interpolation won't work
        TemporalMask(excludeTP)=[];
        WriteInformation(fid,'Excluding first frames from interpolation because there are no non-motion frames before...');
    end
    if any(indd==nTP)
        diffMask=diff(indd);
        excludeTP=indd(find(diffMask>2,1,'last')+1):nTP; % I choose 2 as a threshold here, in case one scan remains between excluded scans, here I will skip this one as well because interpolation won't work
        TemporalMask(excludeTP)=[];
        WriteInformation(fid,'Excluding last frames from interpolation because there are no non-motion frames afterwards...');
    end
    
    
    
    % TCon has size n_vox x n_TP and contains the values of the data points
    % that we know (do not interpolate)
    TCon = fData(:,logical(TemporalMask));

    % tinter is all the time points
    tinter = 1:length(TemporalMask);

    % torig is the time points that we know
    torig = tinter(logical(TemporalMask));

    % We interpolate all time point values (tinter) using the info that we
    % know (torig,TCon)
    % TC has size n_vox x n_TP
    TC = interp1(torig,TCon',tinter,param.interType)';

    WriteInformation(fid,[param.interType ' interpolation finished successfully...']);
end