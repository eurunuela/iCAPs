%% This function constructs a filter given its zeros or its poles; the
% employed equations are derived in 'A Signal Processing Approach to
% Generalized 1-D Total Variation', p. 5267 (discrete versions of
% continuously defined filters)
%
% Inputs:
% - root is a vector containing the zeros or the poles for the filter
%
% Outputs:
% - fil contains the coefficients of the generated filter
%
% Initial version: 20.12.2010, Isik Karahanoglu
function fil = cons_filter(root)

    n=length(root);
    fil = zeros(1,n+1);
    fil(1)=1;
    
    for i = 1:n;
        fil(i+1) = (-1)^i*sum(exp(sum(nchoosek(root,i),2)));
    end

end