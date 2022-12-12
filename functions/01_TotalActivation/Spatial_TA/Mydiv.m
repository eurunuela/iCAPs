% computes a vectorized divergence form GSPbox output
% by Younes
function [out]=Mydiv(s,G)

	tmp = sparse(G.v_in, G.v_out, s, G.N,G.N);
	tmp = tril(tmp,1) - tril(tmp,1)';
	out = sum(tmp.* sqrt(G.W), 2);

end