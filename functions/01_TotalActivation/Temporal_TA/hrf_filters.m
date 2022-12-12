%% This function generates the filters linked to the hemodynamic response
% function of the system: the reconstruction filter converts neural activation to
% BOLD, while the analysis filter converts a BOLD signal into neural
% activity ('SPIKE'), or into the innovation signal ('BLOCK', combination
% of the deconvolution and derivation steps)
%
% Inputs:
% - TR is the TR of the considered data (in s)
% - condition is the type of desired reconstruction filter: solely mapped
% to deconvolution ('SPIKE') or including derivation ('BLOCK')
% - condition2 is the type of hemodynamic response function that we assume:
% either 'bold' (Friston et al., ??, ??), or 'spmhrf' (??? et al., ??, ??)
%
% Outputs:
% - filter_analyze
% - filter_reconstruct
% - maxeig
%
% Implemented by Isik Karahanoglu, 20.12.2010
function [param] = hrf_filters(param)

    % Type of hemodynamic response function to use
    switch param.HRF
        
        %%%%%%%%%%
        % Add details
        case 'bold'

            eps = 0.54;
            ts = 1.54;
            tf = 2.46;
            t0 = 0.98;
            alpha = 0.33;
            E0 = 0.34;
            V0 = 1;
            k1 = 7*E0;
            k2 = 2;
            k3 = 2*E0 - 0.2;

            c = (1 + (1 - E0)*log(1 - E0)/E0)/t0;

            % zeros
            a1 = -1/t0;
            a2 = -1/(alpha*t0);
            a3 = -(1 + 1i*sqrt(4*ts^2/tf - 1))/(2*ts);
            a4 = -(1 - 1i*sqrt(4*ts^2/tf - 1))/(2*ts);

            % pole
            psi = -((k1+k2)*((1 - alpha)/alpha/t0 - c/alpha) - (k3 - k2)/t0)/(-(k1 + k2)*c*t0 - k3 + k2);

        %%%%%%%%%%
        % Add details about that option
        case 'spmhrf'

            %%%%%%%%%%%%%%%%%%
            % CONSTRUCT MODEL
            % FilZeros = [a1,a2,a3,a4]*TR;
            % FilPoles = psi*TR;
            % FilPoles = [];


            % CONSTRUCT MODEL (SPM_HRF)
            % FilZeros = [a1*4,a1*4,real(a3)*0.7+i*imag(a3)*0.55,real(a4)*0.7+i*imag(a4)*0.55]*TR;
            % FilPoles = psi*TR;
            a1 = -0.27;
            a2 = -0.27;
            a3 =-0.4347-1i*0.3497;
            a4 = -0.4347+1i*0.3497;
            psi = -0.1336;
            %%%%%%%%%%%%%%%%%%

        otherwise
            error('Unknown filter'); 
    end

    % Converts the zeros and the poles into [s]
    FilZeros = [a1,a2,a3,a4]*param.TR;
    FilPoles = psi*param.TR;

    % Builds the discrete filters in the time domain according to Karahanoglu
    % et al. 2011 (p.5267). hnum is the filter with the zeros of the linear
    % differential operator, while h_dc and h_dnc are the causal and non causal
    % parts for its poles
    cons=1;
    hnum = cons_filter(FilZeros)*cons;
    hden = cons_filter(FilPoles);

    % Selects the causal and non-causal parts of the operator for the poles
    causal = FilPoles(real(FilPoles)<0);
    n_causal = FilPoles(real(FilPoles)>0);

    % Shortest Filter, 1st order approximation
    h_dc = cons_filter(causal);
    h_dnc = cons_filter(n_causal);

    % Both causal and non-causal parts of the filter are given back as two
    % elements from the h_d array
    h_d{1} = h_dc;
    h_d{2} = h_dnc;

    % Reconstruction filter construction
    param.filter_reconstruct.num = hnum;
    param.filter_reconstruct.den = h_d;

    % One more zero to the filter
    FilZeros2 = [FilZeros,0];

    % Shortest Filter, 1st order approximation
    hnum2 = cons_filter(FilZeros2)*cons;

    % 1024-element frequency response of the filter in d1, and
    % computation of the maximal eigenvalue
    d1 = freqz(hnum2,hden,1024);
    param.MaxEig = max(abs(d1).^2);

    % In the 'block' case, the analysis filter has one more zero
    % compared to the reconstruction filter
    param.f_Analyze.num = hnum2;
    param.f_Analyze.den = h_d;
end



