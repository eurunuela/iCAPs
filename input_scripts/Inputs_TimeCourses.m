% iCAPs-related information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if set to 1, regression will be forced to run, even if already has been done
param.force_Regression=1;

% which type of regression should be done, 
%   'unconstrained' - as in the original paper [Karahanoglu et al., NatComm 2015]
%   'transient-informed' - recommended (default), see [Zoeller et al, IEEE TMI 2018]
param.regType='transient-informed';

% parameter for soft cluster assignment in transient-informed regression,
% can be only one value or a vector of multiple values
% for details see [Zoeller et al., IEEE TMI 2018]
param.softClusterThres=[1:0.2:2];

% choose if for the evaluation of soft cluster assignment factor the
% correlation between measured and estimated amplitudes should be computed
% and evaluated. otherwise, only the BIC and AIC results will be evaluated
param.evalAmplitudeCorrs=0;

% threshold above which a z-scored iCAPs time course will be considered
% "active" - according to Karahanoglu et al, NatComm 2015 and Zoller et
% al., IEEE TMI 2018 we select the default z-score of |1|
param.activityThres=1;