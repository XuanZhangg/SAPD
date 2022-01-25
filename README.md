# Robust-Primal-Dual-Method-for-Computing-Saddle-Point
The code includes the experiment part in the paper Robust-Primal-Dual-Method-for-Computing-Saddle-Point https://arxiv.org/abs/2111.12743 for distributionally logistic regression with l2 regularizer. The details of Lipschitz constants, regularizer, training rate and other parameters can be found in the paper.

The algorithms include: 
1. SAPD: https://arxiv.org/abs/2111.12743
2. OGDA: https://arxiv.org/abs/2002.05683
3. SMP: https://arxiv.org/abs/0809.0815
4. SMD: https://www.semanticscholar.org/paper/Robust-Stochastic-Approximation-Approach-to-Nemirovski-Juditsky/96167ed3ebc9a2c3270f6ae96043e6f086eed4de

Implemention:
The code is implemented by matlab and requires cvx package http://cvxr.com/cvx/.

1. generate_data.m: load and process data, set parameters.
2. main_''algorithm_name''.m: run algorithms and save the result to data folder.
3. plot_tools/main_plot_''figure type''.m: plot the figure and save it to figure folder.
