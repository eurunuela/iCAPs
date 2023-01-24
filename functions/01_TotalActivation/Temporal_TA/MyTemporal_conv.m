%% This function performs the temporal regularization part of total 
% activation, voxel after voxel
%
% Inputs:
% - TCN is the n_time_points x n_ret_voxels 2D matrix of data input to the
% regularization
% - param is a structure containing all TA-relevant parameters; here, we
% will need the fields 'Dimension' (X, Y, Z and T sizes), 'NbrVoxels'
% (number of voxels to consider for regularization), 'LambdaTempCoef' (used
% to compute regularization coefficients)
%
% Outputs:
% - TC_OUT is the n_time_points x n_ret_voxels 2D matrix of outputs from
% the regularization step
% - param is the updated structure with relevant parameters for TA; added
% fields are 'LambdaTempFin' (vector with final regularization estimates
% for each voxel), 'NoiseEstimateFin' (final noise estimate for each voxel)
%
% Implemented by Isik Karahanoglu, 28.11.2011
% CPU / Mex version and conv options by Etienne Orliac (SCITAS) and Younes Farouj (MIP:lab)

function [TC_OUT,param] = MyTemporal(TCN,param)

% 28.11.2011 isik
% computes temporal regularization for all voxels

    TC_OUT  = zeros(param.Dimension(4),param.NbrVoxels);
    param.LambdaTemp = zeros(param.NbrVoxels,1);

    %paramOUT  = cell(1,param.NbrVoxels);    

    %switch lower(param.METHOD_TEMP)

     % case{'spike','block'}

        if (param.COST_SAVE)
            costtemp = zeros(param.NitTemp,param.NbrVoxels);
        end
        
	%EO: Use parfor = 1 to trigger parallel MEX execution
	%        parfor = 0 to default implementation
	%----------------------------------------------------

	if(isempty(gcp('nocreate')) && param.use_parfor)
	    parpool(24)
	end

	Nbr = param.NbrVoxels;

        if (param.use_parfor==1)
	  
	  LambdaTemp    = zeros(param.NbrVoxels,1);
          LambdaTempFin = zeros(param.NbrVoxels,1);

	  %ticBytes(gcp);

	  parfor i=1:Nbr

            paramOUT2 = zeros(2,1);
            TC_OUT2   = zeros(param.Dimension(4),1);


            coef=cconv(TCN(:,i),param.daub);
            

            LambdaTemp(i) = mad(coef,1)*param.LambdaTempCoef;

            param2 = struct;
            param2.f_Analyze     = param.f_Analyze;
            param2.MaxEig        = param.MaxEig;
            param2.Dimension     = param.Dimension;
            param2.NitTemp       = param.NitTemp;
            param2.VxlInd        = i;
            param2.LambdaTemp    = LambdaTemp(i);
            param2.COST_SAVE     = param.COST_SAVE;

	    lambda =  LambdaTemp(i);
            if (isfield(param,'NoiseEstimateFin') && length(param.NoiseEstimateFin) >= i) 
		lambda = param.NoiseEstimateFin(i);
            end

            Temporal_TA_MEX(TCN(:,i),                    ...
			    param.f_Analyze.num,         ...
			    LambdaTemp(i),               ...
			    param.MaxEig,                ...
			    int32(param.Dimension(4)),   ...
			    int32(param.NitTemp),        ...
			    length(param.f_Analyze.den), ...
			    param.f_Analyze.den{1},      ...
			    param.f_Analyze.den{2},      ...
			    lambda,                      ...
			    param.COST_SAVE,             ...
			    TC_OUT2,                     ...
			    paramOUT2);

            LambdaTempFin(i)    = paramOUT2(1);
            NoiseEstimateFin(i) = paramOUT2(2);
            TC_OUT(:,i) = TC_OUT2;

            if (param2.COST_SAVE)
               costtemp(:,i) = paramOUT2.CostTemp;
            end            
         end

         %tocBytes(gcp)
	 
         param.LambdaTemp       = LambdaTemp;
         param.LambdaTempFin    = LambdaTempFin;
         param.NoiseEstimateFin = NoiseEstimateFin;

       else

         for i=1:param.NbrVoxels,
            
            [coef,  len]  = wavedec(TCN(:,i),1,'db3');
            coef(1:len(1)) = [];

            param.LambdaTemp(i) = mad(coef,1)*param.LambdaTempCoef;
            param.VxlInd        = i;

            [TC_OUT(:,i),paramOUT] = TA_Temporal(TCN(:,i),param);

            param.LambdaTempFin(i)    = paramOUT.LambdasTempFin;
            param.NoiseEstimateFin(i) = paramOUT.NoiseEstimateFin;

            if (param.COST_SAVE)
                costtemp(:,i) = paramOUT.CostTemp;
            end
        end

        if (param.COST_SAVE)
            param.cost_TEMP = [param.cost_TEMP; costtemp];
        end
      end

    end



