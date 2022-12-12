function [TA_on_real,TA_on_surrogate,Thresholding_done] = Check_TA_Files(Path,thresholdingTitle)

if nargin < 2
    thresholdingTitle=[];
    Thresholding_done=nan;
end


if exist(fullfile(Path,'TotalActivation','Activity_related.mat'),'file') && ...
        exist(fullfile(Path,'TotalActivation','Activity_inducing.mat'),'file') && ...
        exist(fullfile(Path,'TotalActivation','Innovation.mat'),'file') && ...
        exist(fullfile(Path,'TotalActivation','param.mat'),'file')
    TA_on_real = 1;
else
    TA_on_real=0;
end

if exist(fullfile(Path,'Surrogate','Activity_related_surrogate.mat'),'file') && ...
        exist(fullfile(Path,'Surrogate','Activity_inducing_surrogate.mat'),'file') && ...
        exist(fullfile(Path,'Surrogate','Innovation_surrogate.mat'),'file') && ...
        exist(fullfile(Path,'Surrogate','param.mat'),'file')
    TA_on_surrogate = 1; % Surrogate data exists and TA has been run
elseif exist(fullfile(Path,'Surrogate','Surrogate.mat'),'file')
    TA_on_surrogate = 2; % surrogate data computed, but TA not run yet
else
    TA_on_surrogate=0;
end

if ~isempty(thresholdingTitle)
    if exist(fullfile(Path,'Thresholding',thresholdingTitle,'SignInnov.mat'),'file') && ...
        exist(fullfile(Path,'Thresholding',thresholdingTitle,'param.mat'),'file')
        Thresholding_done = 1;
    else
        Thresholding_done = 0;
    end
end




