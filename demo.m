clear
clc;
load('FGNETmix2.mat')
warning off

r=82;d=1024;n=1002;
U0 = rand(d,r)*aa*2-aa;
V0 = rand(r,n)*aa;
param.maxiter = 100;
param.OriX = X_Ori;
param.InU = U0;
param.InV = V0;
param.k = 6;
param.tol = 1.0e-2;
[label, model, OutU, OutV, llh, converged, t] = BMD_MoG(X_Noi,r,param);