%% This function computes a decision matrix W for TV minimization
% If two neighbour voxels are in different areas with respect to activation 
% (e.g GM  vs WM+CSF) then their finite difference is not considered in TV 
% minimization
% 
% Inputs:
% - Map is the probabilistic map volume
% - sigma determines how sharp the weights are
% - t is the number of time courses
%
% Outputs:
% - wx, wy, wz are the weight volumes (along X, Y and Z)
%
% Implemented by Younes Farouj, 29.04.2016
function [param] = get_weight_gradient(param)

    % dx, dy and dz are 3D matrices depicting the difference between two elements of the GM
    % map along the X, the Y and the Z directions, respectively. Zero-padding
    % is done at the end of the volumes for each case.
    [dx,dy,dz] = gradient3D_full(param.GM_map);

    % If we change from GM to another tissue type, dx is large and so wx_temp
    % becomes low (low weight). Conversely, if we consider two GM voxels
    % neighbouring each other, the weight will be large
    wx_temp=exp(-abs(dx)./param.sigma);
    wy_temp=exp(-abs(dy)./param.sigma);
    wz_temp=exp(-abs(dz)./param.sigma);

    % Repeats the 'gradient volume' t times, with t the number of time points
    param.weight_x = repmat(wx_temp,[1 1 1 param.Dimension(4)]);
    param.weight_y = repmat(wy_temp,[1 1 1 param.Dimension(4)]);
    param.weight_z = repmat(wz_temp,[1 1 1 param.Dimension(4)]);
end