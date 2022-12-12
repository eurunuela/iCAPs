function [Aggregating_done,Clustering_done,ConsensusClustering_done,Regression_ti_done,Regression_unc_done] = Check_iCAPs_Files(data_Path,iCAPs_Path,cons_Path,reg_Path)

if ~isempty(data_Path)
    if exist(fullfile(data_Path,'AI.mat'),'file') && ...
            exist(fullfile(data_Path,'I_sig.mat'),'file') && ... % I have replaced 'Data' by 'I_sig' in the new toolbox version
            exist(fullfile(data_Path,'final_mask.mat'),'file') && ...
            exist(fullfile(data_Path,'subject_labels.mat'),'file') && ...
            exist(fullfile(data_Path,'time_labels.mat'),'file')
        Aggregating_done = 1;
    else
        Aggregating_done = 0;
    end
else
    Aggregating_done=nan;
end

if nargin > 1 && ~isempty(iCAPs_Path)
    if exist(fullfile(iCAPs_Path,'iCAPs.mat'),'file') && ...
            exist(fullfile(iCAPs_Path,'IDX.mat'),'file')
        Clustering_done = 1;
    else
        Clustering_done = 0;
    end
else
    Clustering_done=nan;
end

if nargin > 2 && ~isempty(cons_Path)
    if exist(fullfile(cons_Path,'AUC.mat'),'file') && ...
            exist(fullfile(cons_Path,'CDF.mat'),'file')
        ConsensusClustering_done = 1;
    else
        ConsensusClustering_done = 0;
    end
else
    ConsensusClustering_done=nan;
end


if nargin > 3 && ~isempty(reg_Path)
    if (exist(fullfile(reg_Path,'TC.mat'),'file') && ...
            exist(fullfile(reg_Path,'TC_stats.mat'),'file')) 
        Regression_ti_done = 1;
    else
        Regression_ti_done = 0;
    end
    if (exist(fullfile(reg_Path,'TC_unc_.mat'),'file') && ...
            exist(fullfile(reg_Path,'TC_unc_stats.mat'),'file'))
        Regression_unc_done = 1;
    else
        Regression_unc_done = 0;
    end
else
    Regression_ti_done=nan;
    Regression_unc_done=nan;
end

