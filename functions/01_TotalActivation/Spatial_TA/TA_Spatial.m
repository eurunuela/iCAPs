function x_out = TA_Spatial(y,param,Op,Adj_Op,evaluate_norm)

% This function computes the TV regularization
%
% F(x) = min ||y - x ||^2 + lambda * ||TV{x}||_1
%
% using the 3D extension of the FISTA:

%Beck, A., & Teboulle, M. (2009). A fast iterative shrinkage-thresholding algorithm for linear inverse problems. 
%SIAM journal on imaging sciences, 2(1), 183-202.





%%%%%%---INPUTS---%%%%%%

% y : 3D or 4D volume % If y is 4D, a 3D algorithm is applied to the first 3 compenenets in parallel

% param : Structure containing the following parameters 


%param.NitSpat 
%param.LambdaSpat : The spatial regularization parameter
%param.Lip : Lipshitz constant  ||transpose(gradient)*gradient||_2^2 % in 3D it is given by 2^2+2^2+2^2=12 
%param.tol  : tolerance between iterations for the stopping criteria

%param.VoxelIdx; %[xi,yi,zi] = nonzero indices to convert 2D x to 4D x_vol
%SIZE=param.Dimension; %the size of the 4D data [v1,v2,v3,t]




%%%%%---OUTPUTS---%%%%%%
% x_out : denoised output
%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Younes, 10.03.2016


x_out = zeros(size(y)); % out

%%%% convert to 4D
x_vol = zeros(param.Dimension);
y_vol = zeros(param.Dimension);
% temp = zeros(SIZE);

for i = 1:length(param.VoxelIdx(:,1))
    y_vol(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:) = y(:,i);
end


p=gcp;
fprintf('There are %d workers in pool.\n', p.NumWorkers);


parfor t=1:param.Dimension(4)
x_vol(:,:,:,t) =MyProx(y_vol(:,:,:,t),G,D,Dt,param)
end

for i=1:length(param.VoxelIdx(:,1))
    x_out(:,i) = x_vol(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:);
end

end
