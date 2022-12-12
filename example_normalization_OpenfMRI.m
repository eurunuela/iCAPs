%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION
% this is an example script if TA has been run in subject space and data
% need to be normalized to MNI before continuing with the clustering part
% of the pipeline
%
% In this script, normalization will be applied using the SPM y*
% deformation fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all;

PathData = 'example data';
Subjects = {'sub-10159','sub-10171'};
thresh_title = 'Alpha_5_95_Fraction_0DOT05';


%% set paths
subj_paths=fullfile(PathData,Subjects);
subj_thres_paths=fullfile(subj_paths,'TA_results','exampleToolbox_openfMRI','Thresholding',thresh_title);
subj_TA_paths=fullfile(subj_paths,'TA_results','exampleToolbox_openfMRI','TotalActivation');
subj_input_paths=fullfile(subj_paths,'TA_results','exampleToolbox_openfMRI','inputData');

out_thres=fullfile(subj_paths,'TA_results','exampleToolbox_openfMRI_MNI','Thresholding',thresh_title);
out_TA=fullfile(subj_paths,'TA_results','exampleToolbox_openfMRI_MNI','TotalActivation');

%% define spm deformation batch (subject-independent fields)
spm_jobman('initcfg');
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';



%% run deformation of SignInnov and Activity_inducing
for iS=1:length(subj_paths)
    mkdir(out_thres{iS});mkdir(out_TA{iS});
    % deformation field
    deform_field=dir(fullfile(subj_paths{iS},'anat','Segmented','y_*'));
    deform_field=fullfile(deform_field.folder,deform_field.name);
    
    % apply to SignInnov
    applyTo=dir(fullfile(subj_thres_paths{iS},'SignInnov.nii'));
    applyTo=fullfile(applyTo.folder,applyTo.name);
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deform_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {applyTo};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    spm_jobman('run',matlabbatch);
    movefile(fullfile(subj_thres_paths{iS},'wSignInnov.nii'),fullfile(out_thres{iS},'SignInnov.nii'))
    
    % apply to Activity_inducing
    applyTo=dir(fullfile(subj_TA_paths{iS},'Activity_inducing.nii'));
    applyTo=fullfile(applyTo.folder,applyTo.name);
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deform_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {applyTo};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    spm_jobman('run',matlabbatch);
    movefile(fullfile(subj_TA_paths{iS},'wActivity_inducing.nii'),fullfile(out_TA{iS},'Activity_inducing.nii'))
    
    % apply to input mask
    applyTo=dir(fullfile(subj_input_paths{iS},'mask.nii'));
    applyTo=fullfile(applyTo.folder,applyTo.name);
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deform_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {applyTo};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    spm_jobman('run',matlabbatch);
    movefile(fullfile(subj_input_paths{iS},'wmask.nii'),fullfile(out_TA{iS},'mask.nii'))
    
    % apply to innovations no-nan mask
    applyTo=dir(fullfile(subj_thres_paths{iS},'mask_nonan.nii'));
    applyTo=fullfile(applyTo.folder,applyTo.name);
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {deform_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {applyTo};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    spm_jobman('run',matlabbatch);
    movefile(fullfile(subj_thres_paths{iS},'wmask_nonan.nii'),fullfile(out_thres{iS},'mask_nonan.nii'))
    
end

%% save a new param file for Thresholding and Acitivity inducing with updated information
for iS=1:length(subj_paths)
    
    fHdr=cbiReadNiftiHeader(fullfile(out_thres{iS},'SignInnov.nii'));
    SignInnov_4D=cbiReadNifti(fullfile(out_thres{iS},'SignInnov.nii'));
    mask_nonan_3D=cbiReadNifti(fullfile(out_thres{iS},'mask_nonan.nii'));
    mask_nonan_3D=~isnan(mask_nonan_3D)&mask_nonan_3D~=0;
    mask_3D=cbiReadNifti(fullfile(out_TA{iS},'mask.nii'));
    mask_3D=~isnan(mask_3D)&mask_3D~=0;
    AI_4D=cbiReadNifti(fullfile(out_TA{iS},'Activity_inducing.nii'));
    
    
    %% saving modified param and SignInnov (Thresholding)
    load(fullfile(subj_thres_paths{iS},'param.mat'));
    
    param.mask=reshape(mask_3D,[],1);
    param.Dimension(1)=size(mask_3D,1);
    param.Dimension(2)=size(mask_3D,2);
    param.Dimension(3)=size(mask_3D,3);
    
    param.mask_nonan=reshape(mask_nonan_3D,[],1);
    
    % these are fields that are specific for thresholding in subject space
    % and won't be required further by the pipeline
    param=rmfield(param,{'PC','mask_threshold1'});
    
    % getting 2D innovations in MNI
    SignInnov=reshape(SignInnov_4D,[],size(SignInnov_4D,4));
    SignInnov=SignInnov(param.mask_nonan,:)';
    
    % masking normalized data
    SignInnov_4D(~repmat(mask_nonan_3D,1,1,1,size(SignInnov_4D,4)))=nan;
    
    % saving MNI data
    save(fullfile(out_thres{iS},'param.mat'),'param','-v7.3');
    save(fullfile(out_thres{iS},'SignInnov.mat'),'SignInnov','-v7.3');
    
    hdr=cbiReadNiftiHeader(fullfile(out_thres{iS},'SignInnov.nii'));
    cbiWriteNifti(fullfile(out_thres{iS},'SignInnov.nii'),SignInnov_4D,hdr,'float32');
    
    %% saving modified param (Total Activation)
    load(fullfile(subj_TA_paths{iS},'param.mat'));
    
    param.mask=reshape(mask_3D,[],1);
    param.mask_3D=mask_3D;
    param.Dimension(1)=size(mask_3D,1);
    param.Dimension(2)=size(mask_3D,2);
    param.Dimension(3)=size(mask_3D,3);
    param.IND=find(mask_3D);
    param.VoxelIdx=[];
    [param.VoxelIdx(:,1),param.VoxelIdx(:,2),param.VoxelIdx(:,3)]=ind2sub(size(mask_3D),param.IND);
    param.NbrVoxels=length(param.IND);
    
    % these are fields that are specific for TA in subject space
    % and won't be required further by the pipeline
    param=rmfield(param,{'GM_map','fHeader','weight_x','weight_y','weight_z',...
        'LambdaTemp','LambdaTempFin','NoiseEstimateFin'});
    
    % getting 2D innovations in MNI
    Activity_inducing=reshape(AI_4D,[],size(AI_4D,4));
    Activity_inducing=Activity_inducing(param.mask,:)';
    
    % masking normalized data
    AI_4D(~repmat(mask_3D,1,1,1,size(AI_4D,4)))=nan;
    
    % saving MNI data
    save(fullfile(out_TA{iS},'param.mat'),'param','-v7.3');
    save(fullfile(out_TA{iS},'Activity_inducing.mat'),'Activity_inducing','-v7.3');
    
    hdr=cbiReadNiftiHeader(fullfile(out_TA{iS},'Activity_inducing.nii'));
    cbiWriteNifti(fullfile(out_TA{iS},'Activity_inducing.nii'),AI_4D,hdr,'float32');
end

