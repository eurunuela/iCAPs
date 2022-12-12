function save4Dnii(path,subFolder,fname,data,hdrFilename,mask1D,dim3D)

% Inputs:
%   mask, dim3D - only required if data in 2D

hdr=cbiReadNiftiHeader(hdrFilename);

if nargin < 7
    dim3D=hdr.dim(2:4);
end

if nargin < 6 && length(size(data))==2
    mask1D=ones(size(data));
end


if length(size(data))==2
    if size(data,1)<size(data,2)
        warning(['data probably inverted - number of voxels in mask: ' num2str(size(data,1)) ', number of frames ' num2str(size(data,2)) ', transposing data ...'])
        data=data';
    end
    
    

    if length(mask1D)~=size(data,1) % if data contains only within-mask voxels
        if nnz(mask1D)~=size(data,1)
            disp('There are not the right number of voxels in the mask, please check! skipping nii saving ...');
            return
        end
        data_4D=zeros(length(mask1D),size(data,2));
        data_4D(mask1D,:)=data;
    else
        if hdr.dim(2)*hdr.dim(3)*hdr.dim(4)~=size(data,1)
            error('save4Dnii: wrong dimension in existing header');
        end
        data_4D=data;
    end
    
    data_4D=reshape(data_4D,dim3D(1),dim3D(2),dim3D(3),size(data,2));
    data=data_4D;
    clear data_4D
end


if hdr.dim(2)~=size(data,1) || ...
        hdr.dim(3)~=size(data,2) || ...
        hdr.dim(4)~=size(data,3)
    error('save4Dnii: wrong dimension in existing header');
end

hdr.dim(5)=size(data,4);

outDir=fullfile(path,subFolder);
if ~exist(outDir);mkdir(outDir);end;

outFile=fullfile(outDir,[fname '.nii']);

cbiWriteNifti(outFile,data,hdr,'float32');



