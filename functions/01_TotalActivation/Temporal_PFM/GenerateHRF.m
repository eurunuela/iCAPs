% signalModel = 1 for R2* only and signalModel = 1 for R2* with S0.
% hrf_p currently not used but is there for possible future implementation.

function [X_hrf, X_hrf_norm] = GenerateHRF(dt, TE, nscans, hrf_p, r2only)

temp = [];
tempTE = [];

[hrf_SPM,~] = spm_hrf(dt);
max_hrf = max(hrf_SPM);
hrf_SPM(length(hrf_SPM+1):nscans) = 0;

temp = hrf_SPM;

for j = 1 : nscans-1 % Appends hrf into matrix moving it down by a position on each iteration
    foo = [zeros(j,1); hrf_SPM(1:(end-j))];
    temp = [temp foo];
end

if ~isempty(TE)
    for i = 1 : size(TE,2)
        tempTE = [tempTE; -TE(i)*temp];
    end
else
    tempTE = temp;
end

if r2only
    X_hrf = tempTE;
else
    X_hrf = [tempTE eye(size(tempTE,1),size(tempTE,2))]; % Appends identity matrix for S0
end

X_hrf_norm = X_hrf./max_hrf;
X_hrf = X_hrf*tril(ones(nscans));
X_hrf_norm = X_hrf_norm*tril(ones(nscans));
end
