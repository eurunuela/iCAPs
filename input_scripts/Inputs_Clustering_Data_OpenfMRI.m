
% General data information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path where we have our data stored 
param.PathData = 'example data';

% Links towards the data of all subjects to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of subjects on which to run total activation (must be a cell array
% with all group/subject names). This is where the TA folder will be
% created (or looked for) for each subject
param.Subjects = {'sub-10159','sub-10171'};

% Number of subjects considered
param.n_subjects = length(param.Subjects);

% Title that we wish to give to this specific run of the scripts for saving
% data, or that was used previously for first steps and that we wish to
% build on now
param.title = 'exampleToolbox_openfMRI_MNI';

% name of the iCAPs output for this data
% if only a subset of subjects should be included in the clustering, this
% can be useful to save those different runs in different folders
param.data_title = [param.title '_allSubjects'];


% information about which TA data shoul be used for clustering:
% thresholding information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Alpha-level at which to look for significance of innovation signal frames
% (first element is the percentile of the lower threshold - negative
% innovations - and second element the upper threshold one - positive
% innovations)
param.alpha = [5 95];

% Fraction of voxels from the ones entering total activation for a given 
% subject that should show an innovation at the same time point, so that
% the corresponding frame is retained for iCAPs clustering
param.f_voxels = 5/100;

% Title used to create the folder where thresholding data will be saved
param.thresh_title = ['Alpha_',strrep(num2str(param.alpha(1)),'.','DOT'),'_',...
    strrep(num2str(param.alpha(2)),'.','DOT'),'_Fraction_',...
    strrep(num2str(param.f_voxels),'.','DOT')];


