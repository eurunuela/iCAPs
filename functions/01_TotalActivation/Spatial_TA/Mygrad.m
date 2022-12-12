% computes a vectorized gradient form GSPbox output
% by Younes
function [out]=Mygrad(s,G)
out = (s(G.v_in) - s(G.v_out))';
end
