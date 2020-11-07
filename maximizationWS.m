function [model,k,R] = maximizationWS(X,R,model)
k = size(R,2);

nk = sum(R,1);
mu = zeros(1,k);%fix mu to zero
% mu = bsxfun(@times, X*R, 1./nk);
w = nk/size(R,1);

Sigma = zeros(1,k);
sqrtR = sqrt(R);
for i = 1:k
    Xo = bsxfun(@minus,X,mu(i));
    Xo = bsxfun(@times,Xo,sqrtR(:,i)');
    Sigma(i) = Xo*Xo'/nk(i);
    Sigma(i) = Sigma(i)+(1e-6); % add a prior for numerical stability
end

model.mu = mu;
model.Sigma = Sigma;
model.weight = w;
model.alpha = model.alpha;
model.lambda = model.lambda;