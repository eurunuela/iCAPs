function saveRegionTable(param,iCAPs_z,final_mask)

if ~isfield(param,'regTab_thres') || ~isfield(param,'regTab_codeBook') || ~isfield(param,'regTab_atlasFile')
    warning('attempt to write regions table: no threshold, codeBook or atlas file defined, skipping!');
    return
end

load(param.regTab_codeBook); % loading code book, the param.regTab_codeBook should contain a code book format struct with name "codeBook"

nReg=codeBook.num; % number of regios
nClus=size(iCAPs_z,1); % number of clusters

% reading atlas
atlas_hdr=spm_vol(param.regTab_atlasFile);
atlas_vol=spm_read_vols(atlas_hdr);

% mask header
mask_hdr=spm_vol(fullfile(param.outDir_main,'final_mask.nii'));
mask_vol=spm_read_vols(mask_hdr);

% mapping atlas to iCAPs maps
atlas_vol=mapVTV(atlas_vol,atlas_hdr,mask_hdr);

% vectorize and mask atlas
atlas_vol=atlas_vol(:);
atlas_vol=atlas_vol(final_mask);

% open table file and write header
tab_f=fopen(fullfile(param.outDir_iCAPs,['iCAP_z_regions.txt']),'w');
fprintf(tab_f,'iCAP \t Lobe \t Region \t Percentile \t mean z-score \t voxels \n');
% fprintf(tab_f,'iCAP & Lobe & Region & Percentile & mean z-score & voxels \n \\hline\n');


for iC=1:nClus
    % thresholding iCAPs map
    iCAP_thres=iCAPs_z(iC,:)';
    iCAP_thres(iCAP_thres<param.regTab_thres)=0;
    
    % go through all regions and compute the percentages
    nVox=zeros(nReg,1);nVox_perc=nVox;mean_z=nVox;
    for iR=1:nReg
        nVox(iR,1)=nnz(atlas_vol==codeBook.id(iR)&iCAP_thres);
        nVox_perc(iR,1)=nVox(iR,1)/nnz(atlas_vol==codeBook.id(iR))*100;
        mean_z(iR,1)=mean(iCAP_thres(atlas_vol==codeBook.id(iR)&iCAP_thres));
    end
    
    % sort the regions depending on the percentage
    [nVox_perc_sorted,sortID]=sort(nVox_perc,'descend');
    mean_z_sorted=mean_z(sortID);
    nVox_sorted=nVox(sortID);
    
    
    
    for iR=1:nReg
        if nVox_perc_sorted(iR)>10 % if there is any voxel in the region
            if iR==1
                fprintf(tab_f,'%d',iC);
            end
            fprintf(tab_f,'\t %s \t %s \t %.2f \t %.2f \t %d \n',codeBook.rname{sortID(iR)},...
                codeBook.name{sortID(iR)},nVox_perc_sorted(iR),mean_z_sorted(iR),nVox_sorted(iR));
%             fprintf(tab_f,'& %s & %s & %.2f & %.2f & %d \\\\\n',codeBook.rname{sortID(iR)},...
%                 codeBook.name{sortID(iR)},nVox_perc_sorted(iR),mean_z_sorted(iR),nVox_sorted(iR));
        end
    end
    
    fprintf(tab_f,'\n');
%     fprintf(tab_f,'\\hline\n');
end

