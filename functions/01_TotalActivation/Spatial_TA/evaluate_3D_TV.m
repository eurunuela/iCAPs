%% This function computes 3D l1 norm total variation of a 3D or 4D volume U
%
% Inputs:
% - y : 3D or 4D volume % If y is 4D, a 3D algorithm is applied to the 
% first 3 components in parallel
%
% Outputs:
% - total : total variation norm of y
%
% Implemented by Younes Farouj
function total = evaluate_3D_TV(y)

    [dx, dy, dz] = gradient3D_full(y);

    amplitude = sqrt(abs(dx).^2 + abs(dy).^2 + abs(dz).^2);

    total = sum(amplitude(:));
end
