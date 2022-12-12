% TCS: positive significant innovations and negative significant
% innovations for all voxels (tagged with +1 and -1, respectively)
function [out2d] = check_interconnectedness(data2d,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% check sizes
nVox=size(data2d,2);
nTP=size(data2d,1);
if nVox<nTP
    warning(['data probably inverted - number of voxels: ' num2str(nVox) ', number of frames ' num2str(nTP)])
end


data3d=zeros(param.Dimension(1)*param.Dimension(2)*param.Dimension(3),nTP);

data3d(param.mask_nonan,:)=data2d';
data3d=reshape(data3d,param.Dimension(1),param.Dimension(2),param.Dimension(3),size(data2d,1));

if ~isfield(param,'threshold_interconnectivity') || isempty(param.threshold_interconnectivity)
    WriteInformation(fid,'interconnectivity threshold not specified, defining default of 26');
    param.threshold_interconnectivity=26; % default for 3D data
end

if ~isfield(param,'threshold_minclussize') || isempty(param.threshold_minclussize)
    WriteInformation(fid,'minimum cluster size not specified, defining default 6')
    param.threshold_minclussize=6; % threshold: at least 6 connected clusters
end

for iter_tc=1:size(data3d,4)
    CC = bwconncomp(squeeze(data3d(:,:,:,iter_tc)),param.threshold_interconnectivity);
    CCl = find(cellfun(@length,CC.PixelIdxList)<param.threshold_minclussize);
    out3d = data3d(:,:,:,iter_tc);
    for iter_c = 1:length(CCl)
        out3d(CC.PixelIdxList{CCl(iter_c)})=0;
    end
    data3d(:,:,:,iter_tc) = out3d;
end

out2d=reshape(data3d,[],size(data3d,4));
out2d=out2d(param.mask_nonan,:)';
end


