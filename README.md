We now focus on the factor analysis, on a non-linear (neural net) and linear fashion (tensor decomp), without the local association of each gene (cis-).


The `init_simu.py` simply uses PCA and linear-system solver to init the simulated data, while the `init_real_1.py` (for second layer) and `init_real_2.py` (for first layer) use sparsity solver (LASSO and group LASSO) to get more sparse initialization.



