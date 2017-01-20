**G**enome-wide and **T**ranscriptome-wide **M**ixed-layered **N**eural **N**etwork (GTMNN)

We have the full model here, and the scripts to make three tests happen:

1. direct trans- one-hidden-layer neural net
2. cis- linear regression
3. conditional (on cis-) trans- neural net

This dir contains the following scripts:

1. initialization
2. training
3. analysis (including plot)
4. others

Trans- init: The `init_simu.py` simply uses PCA and linear-system solver to init the simulated data, while the `init_real_1.py` (for second layer) and `init_real_2.py` (for first layer) use sparsity solver (LASSO and group LASSO) to get more sparse initialization.

Cis- init: use LASSO as the average count of cis- SNPs for genes are too large

Conditional trans-: use the residuals from gene profiles taking away cis- effects


