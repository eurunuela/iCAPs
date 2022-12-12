%% This function performs the temporal regularization process on one
% provided time course
%
% Inputs:
% - y - noisy signal, one voxel time course
% - n = ParametersIn.f_Analyze.num;
% - d = ParametersIn.f_Analyze.den;
% - maxeig = Param.MaxEig;
% - N = ParametersIn.Dimension(4);
% - Nit = ParametersIn.NitTemp;
% - lambda = ParametersIn.LambdaTemp(ParametersIn.VxlNbr);
% - ParametersIn.VxlNbr = voxel number
%
% Outputs:
% - x: result of denoising
% - ParametersOut.NoiseEstimateInitial = noise_estimate;
% - ParametersOut.NoiseEstimateFin = nv;
% - ParametersOut.Cost = cost;
% - ParametersOut.LambdasTempFin = Lambda;
%
% Implemented by Isik Karahanoglu, 28.11.2011
function [x,ParametersOut] = TA_Temporal_OneTimeCourse(y,idx_vox,ParametersIn)
    
    % Numerator and denominator of the 'deconv + derivation' operator
    n = ParametersIn.filter_analyze.num;
    d = ParametersIn.filter_analyze.den;

    % Maximal eigenvalue of the operator (used to determine the step size of
    % the algorithm for convergence)
    maxeig = ParametersIn.maxeig;

    % Number of time points
    N = ParametersIn.Dimension(4);

    % Current value of the number of iterations to run temporal regularization
    % for
    Nit = ParametersIn.NitTemp;

    % If we have already went through the whole process, though, lambda has
    % been updated into a gradually more accurate estimate, and se we re-use
    % the final value obtained as a start lambda
    if (isfield(ParametersIn,'NoiseEstimateFin') && (length(ParametersIn.NoiseEstimateFin)>=idx_vox))
            lambda = ParametersIn.NoiseEstimateFin(idx_vox);
    % If we are on the first iteration during which we are running the
    % algorithm, then we take the MAD of wavelet coefficients (computed before
    % the call to Temporal_TA) as 'lambda'
    else
        lambda = ParametersIn.LambdaTemp(idx_vox);
    end

    % Noise_estimate essentially contains the MAD of wavelet coefficients (i.e.
    % the initial value of the regularisation parameter at iteration 1 of the
    % call, when we are into the first forward-backward call
    noise_estimate = ParametersIn.LambdaTemp(idx_vox);

    % nv contains our noise estimate made at each iteration of the algorithm
    % (how far is our solution from the recorded data)
    nv = zeros(Nit,1);

    % Lambda will contain the regularization estimates at each iteration
    Lambda = zeros(Nit,1);

    % Arbitrarily set precision threshold
    precision = noise_estimate/100000;

    % z will contain our dual variable
    z = zeros(N,1);

    % k indexes the iteration we lie in
    k = 1;

    % t is an auxiliary variable used in the weight boosting process (to speed
    % up convergence)
    t = 1;

    % s is the other auxiliary variable used in the boosting process (the
    % 'boosted' version of z)
    s = zeros(N,1);

    % We run the optimisation algorithm for Nit iterations; note that there is
    % no other control for whether we converge or not, only this pre-selected
    % value
    while (k <= Nit)

        % z_l contains the dual estimate from iteration (k)
        z_l = z;

        % Computation of the dual variable at iteration (k+1); it involves the
        % regularisation parameter from the former iteration (lambda), the
        % 'step' of the algorithm (maxeig), and the application of the filter
        % (deconvolution + derivative) defined by n and d to the recorded
        % signals, y, or the (k)th boosted dual estimate, s. The last step
        % (second line) is clipping
        z = 1/(lambda*maxeig)*filter_boundary(n,d,y,'normal') + s - filter_boundary(n,d,filter_boundary(n,d,s,'transpose'),'normal')/maxeig;
        z = max(min(z,1),-1);

        % z now contains the dual estimate from iteration (k+1)

        % t_l stores the (k)th value, and t contains the (k+1)th one. The
        % updates of t and s are part of the 'boosting' scheme that enables to
        % speed up convergence of optimisation algorithms under some conditions
        % met in our case (if I am not wrong, the fact that one of the parts 
        %(of the cost function is smooth and convex)
        t_l = t;
        t = (1+sqrt(1+4*(t^2)))/2;
        s = z + (t_l - 1)/t*(z-z_l);

        % nv measures 'how far we are from the recorded data' at the present 
        % iteration: in fact, it is related to 'y - x(k+1)' with x(k+1) our
        % actual estimate. Because we assume that we have an estimate of true
        % activity and that what is not included in x is the noise, this is our
        % effective noise quantification
        nv(k) = sqrt(sum((lambda*filter_boundary(n,d,z,'transpose')).^2)/N);

        % If the effective noise estimate, at the present iteration, is
        % different from our initial noise estimate, then we want to modify the
        % regularisation parameter accordingly: if we have little 'effective
        % noise', nv(k) is smaller than noise_estimate and lambda increases,
        % because right now, we are including some noise in our solution x
        % whereas we should not, and so we want to sparsify more to get rid of
        % it. Conversely, if we have a very large 'effective noise', lambda
        % will decrease, because we need to get closer to the real data
        if(abs(nv(k) - noise_estimate)>precision);
            lambda = lambda*noise_estimate/nv(k);
        end

        % Lambda stores the updated values of lambda one after the other (the
        % first element is thus the first updated one)
        Lambda(k) = lambda;

        % We move to the next iteration
        k=k+1;
    end

    % When we exit the algorithm, we can use our converged dual variable
    % estimate to get an estimate of the primal variable (our x)
    x = y - (lambda)*filter_boundary(n,d,z,'transpose');

    % ParametersOut.NoiseEstimateIn = noise_estimate;

    % Stores, for the considered voxel, the final noise estimate that was made,
    % and the final related regularisation parameter, so that if we re-enter
    % Temporal_TA, they are directly used as initial estimate
    ParametersOut.NoiseEstimateFin = nv(end);
    ParametersOut.LambdasTempFin = Lambda(end);

end








