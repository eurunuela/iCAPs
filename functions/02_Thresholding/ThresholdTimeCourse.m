%% This function determines the moments of an innovation time course that
% are significant, from comparison to threshold values
%
% Inputs:
% - S is a n_tp-long time course of innovation
% - T is a 2-element vector with lower and upper threshold
%
% Outputs:
% - Out is the returned vector ('1' if significant excursion above 
% upper threshold, '-1' if significant excursion below lower threshold, '0'
% else)
function [Out] = ThresholdTimeCourse(S,T)
    
    % check sizes
    nVox=size(S,2);
    nTP=size(S,1);
    
    if size(T,2)~=nVox
        error('Thresholding time courses: wrong number of threshold pairs, one threshold per voxel required!')
    end
    
    if size(T,1)~=2
        error('Thresholding time courses: wrong number of thresholds, one positive and one negative threshold required!')
    end
    
    
    % Initially filling the output with zeros. If we find a data point
    % lying above the thresholds, we change the related value
    Out = zeros(size(S));
    Out(S <= repmat(T(1,:),nTP,1)) = -1;
    Out(S >= repmat(T(2,:),nTP,1)) = 1;
end