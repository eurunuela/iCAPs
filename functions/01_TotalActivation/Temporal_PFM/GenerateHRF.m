%% This function generates the HRF matrix for the temporal regularization with Paradigm Free Mapping
%
% Inputs:
% - tr: TR of the acqyisition
% - nscans: number of scans of the acquisition
% - block: 1 if you want to use the block model, 0 if you want to use the spike model
%
% Outputs:
% - X_hrf_norm: Normalized HRF matrix for the temporal regularization
%
% Implemented by Eneko Uruñuela, 13.12.2022

function [X_hrf_norm] = GenerateHRF(tr, nscans, block, custom)

temp = [];
tempTE = [];

% If custom is a string, it will be used as the name of the file containing the custom HRF
if ischar(custom)
    hrf = readmatrix(custom);
else
    [hrf,~] = spm_hrf(tr);
end
max_hrf = max(hrf);
hrf(length(hrf+1):nscans) = 0;

temp = hrf;

for j = 1 : nscans-1 % Appends hrf into matrix moving it down by a position on each iteration
    foo = [zeros(j,1); hrf(1:(end-j))];
    temp = [temp foo];
end

X_hrf = temp;

X_hrf_norm = X_hrf./max_hrf;

if block == 1
    X_hrf = X_hrf*tril(ones(nscans));
    X_hrf_norm = X_hrf_norm*tril(ones(nscans));
end
end
