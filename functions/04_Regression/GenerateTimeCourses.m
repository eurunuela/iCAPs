%% This function computes the time courses of activity for all the assessed
% iCAPs.
%
% Inputs:
% - Activity_inducing is a n_final_retained_voxels (after removal of NaN) x
% (subj x time) matrix of concatenated
% subject activity-inducing data
% - iCAPs_rear is a n_iCAPs x n_final_retained_voxels matrix
% - param is the parameter structure
%
% Outputs:
% - T_array is a cell array, with each cell one subject data, of size 
% n_iCAPs x n_timepoints
function [T_array,stats] = GenerateTimeCourses(Activity_inducing,AI_subject_labels,iCAPs,param)
    tic 
    
    if size(Activity_inducing,1)<size(Activity_inducing,2)
        warning('Inverting activity inducing signal for unconstrained regression');
        Activity_inducing=Activity_inducing';
    end
    
    % Number of subjects
    n_sub = param.n_subjects;

    % Number of iCAPs
    n_iCAPs = param.K;

    % Number of time points across all subjects
    n_timepoints = size(Activity_inducing,2);

    % Number of time points for one subject
    for iS=1:n_sub
        n_vol(iS) = nnz(AI_subject_labels==iS);
    end
    
    n_vox = size(Activity_inducing,1);

    % Initial condition for the optimization problem solving (zero for each
    % iCAP)
    x0 = zeros(n_iCAPs,1);

    % Other parameters
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];

    % Bounds within which the estimates are allowed to vary for both the
    % positive part estimation (_pos) and negative part one (_neg)
    lb_pos = zeros(n_iCAPs,1);
    ub_pos = Inf(n_iCAPs,1);
    lb_neg = -Inf(n_iCAPs,1);
    ub_neg = zeros(n_iCAPs,1);

    % Type of optimization algorithm used to solve (set similar to Isik's in
    % her first version)
    options = optimset('Algorithm','interior-point','Display','off');

    % Will contain our estimates
    t_pos = NaN(n_iCAPs,n_timepoints);
    t_neg = NaN(n_iCAPs,n_timepoints);

    % normalize the time series within each subject
%     Activity_inducing_norm=zeros(size(Activity_inducing));
%     for iS=1:n_sub
%         % indices of time series from subject iS
%         vols_iS=AI_subject_labels==iS;
%         % divide by standard deviation
%         Activity_inducing_norm(:,vols_iS)=Activity_inducing(:,vols_iS)./...
%             repmat(sqrt(mean(Activity_inducing(:,vols_iS).^2,2)),1,n_vol);
%         % remove mean across voxels
%         Activity_inducing_norm(:,vols_iS)=Activity_inducing_norm(:,vols_iS) - repmat(mean(Activity_inducing_norm(:,vols_iS),1),n_vox,1);
%     end
    
    
    % For each time point, we perform the estimation
    for t = 1:n_timepoints
        t_pos(:,t) = fmincon(@(time) sum(squeeze((iCAPs'*time - Activity_inducing(:,t)).^2)),x0,A,b,Aeq,beq,lb_pos,ub_pos,nonlcon,options);

        t_neg(:,t) = fmincon(@(time) sum(squeeze((iCAPs'*time - Activity_inducing(:,t)).^2)),x0,A,b,Aeq,beq,lb_neg,ub_neg,nonlcon,options);
    end

    T = t_pos + t_neg;

    % Creation of the array to return
    T_array = {};

    % Loop that retrieves the data for subject s and puts it in a cell array
    for iS = 1:n_sub
       T_array{iS} = T(:,AI_subject_labels==iS);
    end
    
    % compute residuals
    res_pos=iCAPs'*t_pos - Activity_inducing;
    res_neg=iCAPs'*t_neg - Activity_inducing;
    
    % compute model statistics
    for iS=1:n_sub
        vols_iS=AI_subject_labels==iS;
        
        stats.RSS(iS,1)=sum(res_pos(vols_iS).^2+res_neg(vols_iS).^2); % residual sum of squares, sum of positive and negative residuals
        stats.n(iS,1)=n_vox*n_vol(iS)*2; % number of observations (voxels*frames), times 2 (because we do two regressions)
        stats.k(iS,1)=n_iCAPs*n_vol(iS)*2; % number of regressors = number of time points in the unconstrained case, times 2 (because we do two regressions)
        stats.bic(iS,1)=stats.n(iS)*log(stats.RSS(iS)/stats.n(iS)) + stats.k(iS)*log(stats.n(iS));
        stats.aic(iS,1)=stats.n(iS)*log(stats.RSS(iS)/stats.n(iS)) + stats.k(iS)*2;
    end
    
    toc

end
