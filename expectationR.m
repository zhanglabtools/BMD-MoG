function [R, llh] = expectationR(X, model, U, V)
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;
alpha = model.alpha;
lambda = model.lambda;
n = size(X,2);
k = size(mu,2);
logRho = zeros(n,k);

for i = 1:k
    logRho(:,i) = loggausspdf(X,mu(i),Sigma(i));
end
logRho = bsxfun(@plus,logRho,log(w));

T = logsumexp(logRho,2);
E = (diag(alpha-1))*log(V);
F = -sum(sum(abs(U)))/lambda;
llh = sum(T)+sum(sum(E))+F;
llh =llh/n;% loglikelihood
logR = bsxfun(@minus,logRho,T);
R = exp(logR);