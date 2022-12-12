%% This function computes the 3D gradient of a 3D or 4D volume
% It is particulary suitable for fMRI data,  the 4th dimension of V should 
% refer to time
%
% Inputs:
% - V : 3D or 4D volume; if y is 4D, a 3D algorithm is applied to the first 
% 3 components
% - wx, wy, wy : Weight matrices along the 3 variables
%
% Outputs:
% - dx dy dz : gradient
%
% Implemented by Younes Farouj, 12.03.2016
% Weight option added by Younes Farouj, 28.04.2016
function [dx,dy,dz] = gradient3D_full(V,wx,wy,wz)

    % Computes the difference between neighboring elements of the GM map, and
    % pads with a zero at the end of each dimension

    % dx is a 3D matrix (if the GM map is passed as an input argument) with 
    % gradient along X, and same for dy and dz along Y and Z
    dx = padarray( V(2:end,:,:,:) - V(1:end-1,:,:,:), [1 0 0 0], 'post');
    dy = padarray( V(:,2:end,:,:) - V(:,1:end-1,:,:), [0 1 0 0], 'post');
    dz = padarray( V(:,:,2:end,:) - V(:,:,1:end-1,:), [0 0 1 0], 'post');

    % Only used if more than one argument is passed (not the case for now)
    if nargin>1
        dx = dx .* wx;
        dy = dy .* wy;
        dz = dz .* wz;
    end
end 
