function [IDX_mat,dist_thres,clos_to_centroid]=getIDXmat(dist_to_centroid,softClusterThres)

nClus=size(dist_to_centroid,2);

% compute the minimum distance to a cluster for every frame
mindist=min(dist_to_centroid,[],2);
% the threshold distance will be x percent longer than the minimum
% distance (softClusterThres should be a value larger than one)
threshold_dist=softClusterThres*mindist;
threshold_dist=repmat(threshold_dist,1,nClus);

% if a frame is closer to another cluster than the threshlold
% distance, the frame will also be included in this other
% cluster
IDX_mat=zeros(size(threshold_dist));
IDX_mat(dist_to_centroid<=threshold_dist)=1;

dist_thres=dist_to_centroid;
dist_thres(dist_thres>threshold_dist)=0;

clos_to_centroid=1./dist_thres;
clos_to_centroid(clos_to_centroid==Inf)=NaN;
clos_to_centroid=clos_to_centroid./repmat(nansum(clos_to_centroid,2),1,size(clos_to_centroid,2));

