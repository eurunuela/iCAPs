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

function [X_hrf_norm] = GenerateHRF(tr, nscans, block)

temp = [];
tempTE = [];

[hrf_SPM,~] = spm_hrf(tr);
max_hrf = max(hrf_SPM);
hrf_SPM(length(hrf_SPM+1):nscans) = 0;

temp = hrf_SPM;

for j = 1 : nscans-1 % Appends hrf into matrix moving it down by a position on each iteration
    foo = [zeros(j,1); hrf_SPM(1:(end-j))];
    temp = [temp foo];
end

X_hrf = temp;

X_hrf_norm = X_hrf./max_hrf;

if block == 1
    X_hrf = X_hrf*tril(ones(nscans));
    X_hrf_norm = X_hrf_norm*tril(ones(nscans));
end
end
