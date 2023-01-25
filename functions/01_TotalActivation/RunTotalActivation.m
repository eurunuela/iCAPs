%% This function is the core script for total activation, calling the
% spatial and temporal subparts of the routine
%
% Inputs:
% TCN (n_time_points x n_ret_voxels) is the matrix of input data
% param is the structure containing relevant TA parameters; here, we
% require the fields 'Dimension' (4-element vector of X, Y, Z and T sizes),
% 'NbrVoxels' (number of voxels entering TA), 'GM_map' (3D volume
% containing the elements from the probabilistic gray matter map), 'sigma'
% (used in the computation of the spatial TA weights), 'NitTemp' (number of
% iterations to run the temporal scheme for), 'Nit' (number of iterations
% to run the forward-backward scheme for), 'weights' (weight given to
% temporal and spatial regularization schemes in the final averaging
% process)
%
% Outputs:
% - TC_OUT is a n_time_points x n_ret_voxels 2D matrix of data containing
% the output time courses from the total activation process
% Younes, Code optimization, Oct.2019
function [TC_OUT,param] = RunTotalActivation(TCN,param)

    % Spatial regularization is called (note that temporal regularization is
    % called from within MySpatial); TCN contains the time courses that we
    % measure (y), atlas contains the atlas map, param contains all parameters
    % for the algorithm (as well as the gray matter map)
    % Important: if using TV_PPM as a spatial method, 'atlas' is NOT USED
    % anywhere, although it is an argument of MySpatial
    
    TC_OUT = zeros(param.Dimension(4),param.NbrVoxels);

    % Initialization of temporal and spatial matrices (solutions to the
    % two regularization problems)
    xT = zeros(param.Dimension(4),param.NbrVoxels);
    xS = zeros(param.Dimension(4),param.NbrVoxels);

    k=1;
    stepsize=1;

    % The outputs of get_weight_gradient are 4D matrices (time is the
    % fourth dimension, volume the first three): large weights mean
    % that along the x, y or z direction, a given voxel has the same
    % tissue type as his neighbour
    






    
    if ~isfield(param, 'maxit'), param.maxit = 200; end
    if ~isfield(param, 'nu'), param.nu = 1; end
    if ~isfield(param, 'tol'), param.tol = 10e-4; end

	

    disp('Computed weights, entering loop...');
    
    % To solve the problem, the temporal and spatial regularisations
    % are applied one after the other, and then the output from one
    % iteration is the weighted sum of the two outputs; this is known
    % as the 'generalized forward-backward' splitting scheme
    while (k <= param.Nit)

            fprintf('on iteration %d / %d\n', k, param.Nit);
        
        % Increases the number of iterations over which temporal
        % regularization is run at every call (for convergence)
        param.NitTemp = param.NitTemp+100;

        % 1. TEMPORAL REGULARIZATION
        %
        % TC_OUT is only made of zeros at iteration 1, and then is the
        % current solution. TCN contains our data. xT is made of zeros
        % at first iteration and then, contains the temporal
        % regularizer solution
        % Hence, at iteration k=1, we give our data (TCN) to My
        % Temporal. Then, we give (whole estimate - temporal estimate +
        % our data). See Karahanoglu et al. 2013 (NI), Algorithm 1
        % Use the following line if the wavelet toolbox is available.
        %[temp,param] = TA_Temporal(TC_OUT-xT+TCN,param);
        % Otherwise,
        
        if (length(param.Dimension) ~= 4)
            error('SIZE should have 4 dimensions.')
        end



    
        %EO: store last estimate of noise
        noiseFin = zeros(param.NbrVoxels,1,'double');





        TC_IN = TC_OUT-xT+TCN;

            temp = zeros(size(TC_OUT),'double');
            
            tmt = 0;

            % Use synthesis-based temporal regularization with Paradigm Free Mapping (PFM) if
            % param.use_pfm==1, otherwise use the temporal regularization in Total Activation
            if (param.use_pfm==1)
                fprintf('Launching MyTemporal with Paradigm Free Mapping\n');
                tmt = toc;
                [temp,Activity_inducing,Innovation, param] = MyTemporal_pfm(TC_IN, param);
                fprintf('PFM completed in %.5f\n', toc-tmt);
                if (k==param.Nit)
                    param.Activity_inducing=Activity_inducing;
                    param.innovation=Innovation;
                end
            else
                if (param.use_cuda==0)
                    fprintf('Launching MyTemporal with Total Activation\n');
                    tmt = toc;
                    [temp, param] = MyTemporal_conv(TC_IN, param);
                    fprintf('MyTemporal_conv completed in %.5f\n', toc-tmt);
                else
                    fprintf('Launching MyTemporal_MEX with Total Activation\n');
                    tmt = toc;
                    f = parfeval(@MyTemporal_MEX, 2, TC_IN,     ...
                                 param.f_Analyze.num,           ...
                                 length(param.f_Analyze.den),   ...
                                 param.f_Analyze.den{1},        ...
                                 param.f_Analyze.den{2},        ...
                                 param.LambdaTempCoef,          ...
                                 param.MaxEig,                  ...
                                 param.COST_SAVE,               ...
                                 param.NitTemp,                 ...
                                 noiseFin);
                end 
            end
            
            %fprintf('Temporal \n');
            %xT = xT + (temp - TC_OUT); %update temporal, stepsize=1; xT = xT + stepsize*(temp - TC_OUT);
            


        % 2. SPATIAL REGULARIZATION

        % Exactly the same process is done for the spatial
        % regularization; the if condition forces the algorithm to stop
        % with a temporal regularization step (no spatial
        % regularization done at k=5)




       if(k<param.Nit)
                fprintf('Launching Spatial regularization\n');
                tms = toc;
                % calculates for the whole volume
                temp2 = TA_Spatial_Graph(TC_OUT-xS+TCN,param); % calculates for the whole volume
                fprintf('TA_Spatial_Graph completed in %.5f sec\n', toc-tms);
                xS = xS+(temp2-TC_OUT);
            end

            % EO: Wait for parallel execution of MyTemporal
            %     /!\ Update noiseFin with noiseFinOut!
            if (param.use_cuda==1)
                [temp, noiseFin] = fetchOutputs(f);
                fprintf('MyTemporal_MEX completed in %.5f\n', toc-tmt);
                %noiseFin = noiseFinOut
            end


        % 3. WEIGHTED AVERAGE OF THE SOLUTIONS
        xT = xT + stepsize*(temp - TC_OUT); %update temporal, stepsize=1; xT = xT + stepsize*(temp - TC_OUT);
            
            TC_OUT = xT*param.weights(1)+param.weights(2)*xS;
        disp('Finished weighted averaging step...');

        k = k+1;
    end 
end
