%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DeconvolutionAlgorithm:
%   Extracts neural activity signal from noisy fMRI signal by applying a
%   deconvolution with a FISTA.
%
%   Inputs:
%   - input_signal: struct with data from the signal to be analysed.
%   - r2only: 0 for R2* only, 1 for model with S0.
%   - fista_params: struct with parameters needed for the FISTA.
%       - rho: rho value to weight L1 and L21 in the proximal method.
%       - lambda: lambda values of the proximal method.
%       - weights: weight of lambda values.
%       - itermax: maximum number of FISTA iterations.
%       - converged: convergence criteria to stop FISTA.
%       - proximal_mode: name of the proximal method to be used.
%   - stability: 1 for stability selection purposes, 0 otherwise.
%
%   Outputs:
%   - betastruct: struct with deconvolved activity signal.
%       - beta: non debiased deconvolved activity signal.
%       - debiased: debiased deconvolved activity signal.
%   - Ystruct: struct with bold signal of deconvolved activity signal.
%       - fit: non debiased bold signal.
%       - debiased: debiased bold signal.
%
%   BCBL, July 2018
%   Eneko Urunuela
%   e.urunuela@bcbl.eu

function [betastruct, Ystruct, s_zero] = DeconvolutionAlgorithm(input_signal, fista_params, hrf)

%% Check input parameters or set default values
if (~isfield(fista_params,'proximal_mode')) || isempty(fista_params.proximal_mode)
    fista_params.proximal_mode = 'lasso mixed space';
end
if (~isfield(fista_params,'rho')) || isempty(fista_params.rho)
    fista_params.rho = 0.5;
end
if (~isfield(fista_params,'lambda')) || isempty(fista_params.lambda)
    fista_params.lambda = 100;
end
if (~isfield(fista_params,'itermax')) || isempty(fista_params.itermax)
    fista_params.itermax = 5000;
end
if (~isfield(fista_params,'converged')) || isempty(fista_params.converged)
    fista_params.converged = 1e-12;
end
if (~isfield(fista_params,'weights')) || isempty(fista_params.weights)
    fista_params.weights = [];
end

%% Parameterization
signal          = input_signal.data;
TE              = input_signal.te;
%TR              = input_signal.tr;
TE_mean         = mean(TE);
% TE_norm         = TE/TE_mean; % Not dividing by the mean drastically reduces beta amplitude, reducing the contrast in the imagesc.
% n_hrf           = 1;
nscans          = size(signal,1)/length(TE);
nvoxels         = size(signal,2);
rho             = fista_params.rho;
lambda          = fista_params.lambda;
weights         = fista_params.weights;
itermax         = fista_params.itermax;
converged       = fista_params.converged;
proximal_mode   = fista_params.proximal_mode;

% HRF matrix
X_hrf = hrf.hrf;
X_hrf_norm = hrf.norm;

clear temp

%% Compute variables used in the iterative algorithm
X_tilde = X_hrf_norm;
X_tilde_trans = X_tilde';
X_tilde_tt = X_tilde_trans*X_tilde;
c_ist = 1/abs(eigs(X_tilde_tt,1));
data_tilde = signal;

%% Initialization of variables
v = X_tilde_trans*data_tilde;

beta = zeros(size(X_hrf,2),nvoxels);
Y_fit = zeros(size(signal,1),nvoxels);

% FISTA parameters
y_fista = zeros(size(X_hrf,2),nvoxels);
prox_z = zeros(size(X_hrf,2),nvoxels);

%% FISTA
t_fista     = 1; % Initial value as stated in Gramfort et al 2011

tic

fprintf('Starting FISTA... \n');

for iter = 1 : itermax
    beta_old = beta;
    prox_z_old = prox_z;
    
    y_ista = y_fista;
    
    % Forward-Backward step
    z_ista = y_ista + c_ist*(v - X_tilde_tt*y_ista);
    
    z_ista_hrf = z_ista(1:nscans,:);

    prox_z_hrf = proximal_operator_lasso(z_ista_hrf, c_ist*lambda, weights);
    
    if r2only
        prox_z = prox_z_hrf;
        s_zero = [];
    else
        prox_z = [prox_z_hrf; z_ista(nscans+1:end,:)];
        s_zero = z_ista(nscans+1:end,:);
    end
    
    % Accelerating procedure for ISTA
    beta = prox_z; % Estimates (R2* and S0, or only R2*)
    t_fista_old = t_fista;
    t_fista = 0.5*(1+sqrt(1+4*t_fista_old^2));
    y_fista = prox_z + (prox_z - prox_z_old)*(t_fista_old-1)/t_fista;
    
    % Save estimated bold signals for current iteration
    Y_fit = X_hrf_norm*beta;
    
    % Calculates error between current and previous iteration.
    convergence_criterium = sum(sum((beta - beta_old).^2))/sum(sum(beta_old.^2));
    
    if convergence_criterium <= converged && stability == 0
        fprintf('FISTA has converged! \n');
        break;
    end
    
    fprintf('Iteration %i/%i \n', iter, itermax);

end % FISTA

elapsed_time = toc;

fprintf('Total FISTA time was %.2f seconds \n',round(elapsed_time,2));

% Saves function output as struct
betastruct.beta         = beta;
Ystruct.fit             = X_hrf*beta;
betastruct.debiased     = beta_debiased;
Ystruct.debiased        = X_hrf*kappa;

end % function