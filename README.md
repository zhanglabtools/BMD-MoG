# BMD-MoG
The code file for paper "Robust Bayesian Matrix Decomposition with Mixture of Gaussian Noise". The main function is `BMD_MoG`
``` matlab
function [label,model,OutU,OutV,llh,converged,t] = BMD_MoG(InX,r,param)
% function [label,model,OutU,OutV,llh,converged,t] = BMD_MoG(InX,r,param)
% Perform EM algorithm for fitting the BMD_MoG model.
% Step 1: Max parameters;
% Step 2: Expectation;
% Step 3: Max weighted L2 MF;
% Step 4: Expectation.
%Input:
%   InX: d x n input data matrix
%   r:   the rank
%   param.maxiter: maximal iteration number
%   param.OriX: ground truth matrix
%   param.InU,InV: Initialized factorized matrices
%   param.k: the number of GMM
%   param.display: display the iterative process
%   param.tol: the thresholding for stop
%Output:
%   label:the labels of the noises
%   model:model.mu, the means of the different Gaussians
%         model.Sigma,the variance of the different Gaussians,sigma^2
%         model.weight,the mixing coefficients
%         model.alpha,the parameter of  Dirichlet prior on V
%         model.lambda,the parameter of Laplace prior on U
%   W: d x n weighted matrix
%   OutU: the fianl factorized matrix U
%   OutV: the fianl factorized matrix V
%   llh:  the log likelihood
%
```
