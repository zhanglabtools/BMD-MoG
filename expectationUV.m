function [U, V] = expectationUV(X,R,U,V,model)
k = size(R,2);
[d,n]=size(X);
alpha = model.alpha;
Sigma = model.Sigma;
lambda = model.lambda;

Rnew = zeros(size(R,1),1);
for i=1:k
    Rnew = Rnew+R(:,i)./Sigma(i);
end
B = reshape(Rnew,d,n);

for i=1:d
    U(i,:) = softthresholding(X(i,:)*diag(B(i,:))*V',1/lambda)/(V*diag(B(i,:))*V');
end

for j=1:n
    v0=V(:,j);
    Bj=diag(B(:,j));
    xj=X(:,j);
    alpha_new=alpha-1;
    alpha_new=alpha_new';
   
    Q = U' * Bj * U;
    b = U' * Bj * xj;
        
    cal_objval(alpha_new, b, v0, Q);
    v_new = updateVj(alpha_new, b, v0, Q);
    cal_objval(alpha_new, b, v_new, Q);
    V(:,j)=v_new;
end
end
function [ soft_thresh ] = softthresholding( a,K )
    soft_thresh = sign(a).*max(abs(a) - K,0);
end