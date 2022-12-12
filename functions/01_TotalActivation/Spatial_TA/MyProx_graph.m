% F(x) = min ||y - x ||^2 + lambda * ||TV{x}||_1, where TV is the gradient
% operator using a graph extension of the FISTA:
% Beck, A., & Teboulle, M. (2009). A fast iterative shrinkage-thresholding 
% algorithm for linear inverse problems. 
% SIAM journal on imaging sciences, 2(1), 183-202.
%
% Inputs:
% - y : signal over graph
% - param: structure containing required TA parameters
%
% Outputs:
% - x : denoised output
%
% Implemented by Younes Farouj, 12.02.2018
function [x,nv,J]=MyProx_graph(y,param)




D = param.A;
Dt=param.At;

Nit = param.NitSpat;
lambda=param.LambdaSpat;
stop_cri = param.tol;
maxeig = 2*param.lmax*param.nu;

 


       
N = length(y);
k = 1;
t = 1;
z = zeros(size(D(y)));
s = zeros(size(D(y)));


J = zeros(Nit,1);
nv = zeros(Nit,1);

    
    
  %% main loop  
  
     while (k <=Nit)

        z_l = z;
        z = 1/(lambda*maxeig)*D(y) + s - D(Dt(s)')/maxeig;
        % clipping
        z = max(min(z,1),-1);             
        t_l = t;
        t = (1+sqrt(1+4*(t^2)))/2;
        s = z + (t_l - 1)/t*(z-z_l);  
        x = y - (lambda)*Dt(z)'; % update x
        J(k) = sum(abs(x-y).^2) + lambda*sum(abs(D(x)));
        nv(k) = sqrt(sum((x-y).^2)/N);

%      if(p_fig==1),figure(101);plot(y,'-g');hold on;plot(x,'-r');pause(0.2);hold off;end
         if((k > 3) && (abs(J(k-1)-J(k-2)) < stop_cri)) % stopping criterion 
             k
             break;
         end

         k=k+1;

    end
     x = y - (lambda)*Dt(z)'; % update x
     
end


     
      
      