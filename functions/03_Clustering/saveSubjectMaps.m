function saveSubjectMaps(param,subject_labels,IDX,I_sig,final_mask)
nClus=max(IDX);
nSub=max(subject_labels);
if ~exist(fullfile(param.outDir_iCAPs,'subjectMaps'),'dir');mkdir(fullfile(param.outDir_iCAPs,'subjectMaps'));end;

% create a 4D data matrix with all frames for all subjects
for iC=1:nClus
    iCAP_sub=nan(nSub,size(I_sig,2));
    for iS=1:nSub
        iCAP_sub(iS,:)=mean(I_sig(IDX==iC&subject_labels==iS,:));
    end
    iCAP_sub = ZScore_iCAPs(iCAP_sub);
    save4Dnii(param.outDir_iCAPs,'subjectMaps',['iCAP_z_' num2str(iC)],iCAP_sub',fullfile(param.outDir_main,'final_mask.nii'),final_mask);
end