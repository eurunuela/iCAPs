%% This function computes the following regularization:
% F(x) = min ||y - x ||^2 + lambda * ||Op{x}||_1, where Op is the gradient
% operator using the 3D extension of the FISTA:
% Beck, A., & Teboulle, M. (2009). A fast iterative shrinkage-thresholding 
% algorithm for linear inverse problems. 
% SIAM journal on imaging sciences, 2(1), 183-202.
%
% Inputs:
% - y : 3D or 4D volume % If y is 4D, a 3D algorithm is applied to the 
% first 3 components in parallel
% - param: structure containing required TA parameters
%
% Outputs:
% - x_out : denoised output
%
% Implemented by Younes Farouj, 10.03.2016
function x_out = MyProx(x,G,D,Dt,param)


% Initializing solution and parameters

  nu      = G.lmax*param.nu;
  gamma   = param.gamma;
  stopcri = param.stopcri;
  maxiter = param.maxiter;

 % Initializing solution

  x_out=zeros(size(x));


% Initializing algorithm
  t = 1;
  z = zeros(size(D(x)));
  s = zeros(size(D(x)));

    %% main loop
     for i = 1:maxiter 
        
        x_old = x_out;
            
        z_l = z;
        z = 1/(gamma*nu)*D(x) + s - D(Dt(s))/nu;      % Backward step
        z = max(min(z,1),-1);                         % clipping       
        t_l = t;
        t = (1+sqrt(1+4*t*t)) * 0.5;                  % FISTA Acceleration 
        s = z + (t_l - 1)/t*(z-z_l);  
        x_out = x - gamma*Dt(z);                      % update x         
        
         error = norm(x_out(:) - x_old(:),2)/norm(x_out(:),2);
        if error < stopcri
            break;
        end
        
     end
 
    x_out = x - gamma*Dt(z); % update x

 
end