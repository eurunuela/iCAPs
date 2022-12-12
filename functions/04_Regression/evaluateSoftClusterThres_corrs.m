
function evaluateSoftClusterThres_corrs(clusteringResults,TC,param,fid)

thresVals=param.softClusterThres;
AI=clusteringResults.AI';
AI_subject_labels=clusteringResults.AI_subject_labels;
subject_labels=clusteringResults.subject_labels;
time_labels=clusteringResults.time_labels;
dist_to_centroid=clusteringResults.dist_to_centroid;
iCAPs_z=ZScore_iCAPs(clusteringResults.iCAPs);

nThres=length(thresVals);
nSub=param.n_subjects;
nVox=size(AI,1);
nTP_all=size(AI,2);
nI=size(subject_labels,1);

%% correlation between measured and estimated TC amplitudes
for iT=1:nThres
    [IDX_mat,dist_thres,clos_to_centroid]=getIDXmat(dist_to_centroid,thresVals(iT));
    
    dist_mat_all(:,:,iT)=dist_thres;
    IDX_mat_all(:,:,iT)=IDX_mat;
    nClus_mat_all(:,iT)=sum(IDX_mat,2);
    clus_weights(:,:,iT)=clos_to_centroid;
end


disp('Computing normalized innovations...')
Activity_inducing_norm=zeros(size(AI));
innovations_norm=zeros(size(AI));%zeros(size(AI,1),size(AI,2)-nSub);
for iS=1:nSub
    % indices of time series from subject iS
    vols_iS=AI_subject_labels==iS;
    nTP_sub(iS)=nnz(AI_subject_labels==iS);

    % divide by standard deviation
    Activity_inducing_norm(:,vols_iS)=AI(:,vols_iS)./repmat(sqrt(mean(AI(:,vols_iS).^2,2)),[1,nTP_sub(iS)]);

    % remove mean across voxels (global signal regression)
%         Activity_inducing_norm(:,vols_iS)=Activity_inducing_norm(:,vols_iS) - repmat(mean(Activity_inducing_norm(:,vols_iS),1),n_vox,1);
    Activity_inducing_norm(isnan(Activity_inducing_norm))=0;
    
    innovations_norm(:,vols_iS)=[zeros(size(AI,1),1),diff(Activity_inducing_norm(:,vols_iS),1,2)];   
end


% compute the innovation amplitudes in the clustering maps, threshold 1.5, 0.5 and 2.3
disp('Computing innovation amplitudes in the cluster maps...')
thres=1.5;
disp(['Spatial threshold (zscore): ' num2str(thres)]);

for iT=1:nThres
    disp(['Soft cluster assignment factor: ' num2str(thresVals(iT))]);
    for iI=1:nI
        if ~mod(iI,1000)
            fprintf([num2str(iI) ' ']);
        end
        iS=subject_labels(iI);
        t=time_labels(iI);if t>nTP_sub(iS);t=t-nTP_sub(iS);end
        clusterID=find(IDX_mat_all(iI,:,iT));
        clusterWeights=clus_weights(iI,:,iT);
        vols_iS=AI_subject_labels==iS;

        % compute "innovation amplitude", i.e. the average innovation inside 
        % the mask of the iCAP this frame was assigned to
        innovSub=innovations_norm(:,vols_iS);
        
        if t==1;t=t+1;end
        TC_change_tmp=TC{iT}{iS}(:,t)-TC{iT}{iS}(:,t-1);

        % loop through all clusters
        AI_change_measured(iI,iT)=0;
        TC_change_estimated(iI,iT)=0;
        for iC=clusterID
            AI_change_measured(iI,iT)=AI_change_measured(iI,iT)+...
                    mean(innovSub(innovSub(:,t)~=0&iCAPs_z(iC,:)'>thres,t))*clusterWeights(iC);
            TC_change_estimated(iI,iT)=TC_change_estimated(iI,iT)+TC_change_tmp(iC)*clusterWeights(iC);
        end
    end
end
AI_change_measured(isnan(AI_change_measured))=0;
TC_change_estimated(isnan(TC_change_estimated))=0;

corr_amps=diag(corr(AI_change_measured,TC_change_estimated));

[~,opt_id]=findpeaks(corr_amps);
opt_id=opt_id(1);

WriteInformation(fid,['The optimum correlation is at xi=' num2str(thresVals(opt_id))]);

figure('position',[440   565   270   233]);
plot(thresVals,corr_amps,'*');
hold on
plot(thresVals(opt_id),corr_amps(opt_id),'*r')
title('correlation measured vs. estimated transient amplitudes')
print(fullfile(param.outDir_reg,'evalCorr'),'-depsc2','-painters');










