%% This function regresses out low-frequency components, and possible
% additional covariates, from the data. Low-frequency components are
% obtained from a DCT basis
%
% Inputs:
% - in is the input data (1D time course, n_tp long)
% - TR is the data TR
% - TS is the sampling rate ?
% - covariates contains possible covariates ([] if unwanted)
%
% Outputs:
% dct_sol is the 1D n_tp long solution
% c_dct are the DCT coefficients
function [dct_sol, c_dct] = sol_dct(in,TR,TS,covariates)

    n = length(in);

    k = round(2*n*TR/TS + 1);

    dct = spm_dctmtx(n,k); %construct DCT matrix
    lt = (1:n);
    lt = lt./sqrt(sum(lt.^2)); % linear trend
    %dct = [dct , lt' , covariates];

    c_dct = (transpose(dct)*dct)\transpose(dct)*in;  % DCT coefficients

    dct_sol = in - dct*c_dct; % detrended

end
