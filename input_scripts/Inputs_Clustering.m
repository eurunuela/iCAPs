% iCAPs-related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify if clustering should be done (one may want to only run consensus 
% clustering, then set this to 0), default = 1
param.doClustering=1;

% if set to 1, Clustering will be forced to run, even if already has been done
param.force_Aggregating=0;
param.force_Clustering=0;

% file with external mask to use to define the input voxels,
% if nothing is specified, the mask will be the intersection of the GM
% masks of all subjects that are included in the clustering
param.common_mask_file=[];

% file with additional mask to apply additionally to the intersection of
% all GM masks, it has to be in the same space as all subject's significant
% innovations (I am using this here to exclude the cerebellum)
param.extra_mask_file='GM_mask_MNI333_AAL.nii';




% Number of iCAPs into which to separate the data
param.K = 19:21;

% Type of distance to use for the k-means clustering process (choose
% between 'sqeuclidean' and 'cosine')
param.DistType = 'cosine';

% Number of times the clustering process is run in a row to extract iCAPs
param.n_folds = 10;

% specify if the result of each replicate should be saved during clustering
param.saveClusterReplicateData=0;

% Maximum number of allowed iterations of the kmeans clustering, the Matlab
% default of 100 is sometimes not enough if many frames are included,
% default = 100
param.MaxIter=300;


% save subject-specific iCAPs maps
param.saveSubjectMaps=1;

% save iCAPs regions tables
param.saveRegionTables=1;
param.regTab_thres=1.5; % z-score at which to threshold map for regions table
param.regTab_codeBook='AALcodeBook.mat'; % file with atlas region names
param.regTab_atlasFile='AAL90_correctLR.nii'; % file with atlas data in MNI (has to be in same space as iCAPs results)


% Title used to create the folder where iCAPs data will be saved
for nK=1:length(param.K)
    param.iCAPs_title{nK} = ['K_',num2str(param.K(nK)),'_Dist_',...
            param.DistType,'_Folds_',num2str(param.n_folds)];
end
if length(param.K)==1
    param.iCAPs_title=param.iCAPs_title{1};
end


% consensus clustering - related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select if consensus clustering should be run or not
param.doConsensusClustering=1;

% if set to 1, Consensus Clustering will be forced to run, even if already has been done
param.force_ConsensusClustering=0;

% Subsample Type:
% 'subjects' to subsample all frames from a subject; 
% 'items' to subsample frames without taking into account within- or
%   between-subject information
param.Subsample_type='items';
param.Subsample_fraction=0.8;
param.cons_n_folds=20;
param.cons_title=[num2str(param.K(1)) 'to' num2str(param.K(end)) ...
    '_SubsampleType_' param.Subsample_type ...
    '_Fraction_' strrep(num2str(param.Subsample_fraction),'.','DOT') ...
    '_nFolds_' num2str(param.cons_n_folds) ...
    '_Dist_' param.DistType];


% flag to indicate that cluster consensus should be computed, clustering
% and consensus clustering have to be done to compute this measure
param.computeClusterStability=1;

