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

function [param] = get_decision_gradient(param)

    % dx, dy and dz are 3D matrices depicting the difference between two elements of the GM
    % map along the X, the Y and the Z directions, respectively. Zero-padding
    % is done at the end of the volumes for each case.

[dx,dy,dz] = gradient3D_full(param.GM_map);

wx_temp=1-abs(dx);
wy_temp=1-abs(dy);
wz_temp=1-abs(dz);

    param.weight_x = wx_temp;
    param.weight_y = wy_temp;
    param.weight_z = wz_temp;


end
