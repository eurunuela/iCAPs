function param = make_graph(param)


    G=grid_gm(param);
    G = gsp_adj2vec(G);
    G = gsp_estimate_lmax(G);
    
    param.A  =  @(s) Mygrad(s,G);
    param.At = @(x) Mydiv(x,G);
    param.lmax = G.lmax;


end