%% This function adds all the useful paths to the matlab environment for
% the application of total activation routines
function [] = AddPaths()

    % Path towards the toolbox functions
    addpath(fullfile('functions','00_Utilities'));
    addpath(fullfile('functions','00_Preprocessing'));
    addpath(genpath(fullfile('functions','01_TotalActivation')));
    addpath(fullfile('functions','02_Thresholding'));
    addpath(fullfile('functions','03_Clustering'));
    addpath(fullfile('functions','04_Regression'));
    addpath(fullfile('functions'));
%     addpath(genpath(fullfile('functions')));

    addpath(fullfile('input_scripts'));

    
    % Path towards SPM12
    addpath(fullfile('packages','spm12'));
    
    % Path towards cbiNifti to read and write 4D nifti data
    addpath(fullfile('packages','cbiNifti'));
    
    % Path towards munkres function for Hungarian algorithm
    addpath(fullfile('packages','munkres'));
    
    % Path towards function to find knee points
    addpath(fullfile('packages','knee_pt'));

     % Path towards function to findgsp toolbox
    addpath(genpath(fullfile('packages','gspbox')));
    
end