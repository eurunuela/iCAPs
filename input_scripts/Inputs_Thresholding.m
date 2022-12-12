% Thresholding-related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if set to 1, Thresholding will be forced to run, even if already has been done
param.force_Thresholding=1;


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

% Number of neighbours that must also show an innovation for a voxel to be
% retained
param.threshold_minclussize = 6;

% Number of neighbours to consider in the process
param.threshold_interconnectivity = 26;


