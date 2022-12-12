
% Total activation-related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select if you want to force TA to run, and overwriting existing files if
% it has already been done
param.force_TA_on_real=1;
param.force_TA_on_surrogate=1;


% Type of assumed hemodynamic response function type (can be 'bold' or
% 'spmhrf')
param.HRF ='bold';

% Number of outer iterations for the whole forward-backward scheme:
% this means that at most, the solutions with both regularizer 
% terms will be found and averaged five times. This is the 'k_max'
% in Farouj et al. 2016 (ISBI abstract from Younes)
param.Nit=5;

% Number of iterations for which the temporal regularization scheme
% and the spatial one are run
param.NitTemp = 500;
param.NitSpat=400;

% Tolerance threshold for convergence of the TA methods
param.tol = 1e-5;

% Weighting parameters for the temporal and the spatial TA schemes: because
% the temporal part is more elaborate and more extensively tested, Dr
% Karahanoglu and Dr Farouj opted for a 3/4 to 1/4 trade-off =P
param.weights = [0.75 0.25];

% Coefficient somehow used to multiply the regularization weights of
% temporal regularization
param.LambdaTempCoef = 1/0.8095;


% The noise estimation procedure uses a single scale
% wavelet decomposition. Here Daubechies wavelets with 4 vanishing 
% moments are used. The corresponding high pass filter is given by:

g=[0    -0.12941    -0.22414     0.83652    -0.48296];
g = g';
param.daub = g;


% Weight of spatial regularization
param.LambdaSpat=6;




param.COST_SAVE = 0; % save the costs..


