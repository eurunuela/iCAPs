
function best_id_bic=evaluateSoftClusterThres(TC_stats,param,fid)

thresVals=param.softClusterThres;

nThres=length(thresVals);
nSub=param.n_subjects;


%% BIC/AIC
% get BIC values from stats
BICvals=zeros(nThres,nSub);
AICvals=zeros(nThres,nSub);
for iT=1:nThres
    BICvals(iT,:)=TC_stats{iT}.bic;
    AICvals(iT,:)=TC_stats{iT}.aic;
end

% BIC knee point
[~, knee_idx] = knee_pt(sum(BICvals,2),thresVals);
WriteInformation(fid,['The BIC knee point is at xi=' num2str(thresVals(knee_idx))]);
best_id_bic=knee_idx;

% plotting BIC values
figure('position',[440   560   277   238]);
plot(thresVals,sum(BICvals,2),'*')
hold on;
plot(thresVals(knee_idx),sum(BICvals(knee_idx,:),2),'*r');
title('Sum of BIC across subjects');
print(fullfile(param.outDir_reg,'BICsum'),'-depsc2','-painters');


% plotting BIC values
figure('position',[440   560   277   238]);
errorbar(thresVals,mean(BICvals,2),std(BICvals,[],2))
hold on;
plot(thresVals(knee_idx),mean(BICvals(knee_idx,:),2),'*r');
title('BIC distribution across subjects');
print(fullfile(param.outDir_reg,'BICdist'),'-depsc2','-painters');

% BIC knee point
[~, knee_idx] = knee_pt(sum(AICvals,2),thresVals);
WriteInformation(fid,['The AIC knee point is at xi=' num2str(thresVals(knee_idx))]);

% plotting BIC values
figure('position',[440   560   277   238]);
plot(thresVals,sum(AICvals,2),'*')
hold on;
plot(thresVals(knee_idx),sum(AICvals(knee_idx,:),2),'*r');
title('Sum of AIC across subjects');
print(fullfile(param.outDir_reg,'AICsum'),'-depsc2','-painters');


% plotting BIC values
figure('position',[440   560   277   238]);
errorbar(thresVals,mean(AICvals,2),std(AICvals,[],2))
hold on;
plot(thresVals(knee_idx),mean(AICvals(knee_idx,:),2),'*r');
title('AIC distribution across subjects');
print(fullfile(param.outDir_reg,'AICdist'),'-depsc2','-painters');









