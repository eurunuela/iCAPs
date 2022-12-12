%% This function computes the 3D divergence of the vector field [dx,dy,dz]
% It is particulary suitable for fMRI data,  the 4th dimension of d* should 
% refer to time
%
% Inputs:
% - [dx,dy,dz] : 3D vector field, typically a gradient. If d* is 4D, a 3D 
% algorithm is applied to the first 3 components
% - wx, wy, wy : Weight matrices along the 3 variables
%
% Outputs:
% - V : divergence
%
% Implemented by Younes Farouj, 12.03.2016
% Weight option added by Younes Farouj, 28.04.2016
function V = div3D_full(dx,dy,dz,wx,wy,wz)

    if nargin > 3
        dx = dx .* wx;
        dy = dy .* wy;
        dz = dz .* wz;
    end

    V = padarray( dx(1:end-1,:,:,:), [1 0 0 0], 'post' ) - ...
        padarray( dx(1:end-1,:,:,:), [1 0 0 0], 'pre' ) + ...
        padarray( dy(:,1:end-1,:,:), [0 1 0 0], 'post' ) - ...
        padarray( dy(:,1:end-1,:,:), [0 1 0 0], 'pre' )  + ...
        padarray( dz(:,:,1:end-1,:), [0 0 1 0], 'post' ) - ...
        padarray( dz(:,:,1:end-1,:), [0 0 1 0], 'pre' );

end 
