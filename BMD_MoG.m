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


%% initialization
[d,n] = size(InX);
if (~isfield(param,'maxiter'))
    maxiter = 100;
else
    maxiter = param.maxiter;
end

if (~isfield(param,'OriX'))
    OriX = InX;
else
    OriX = param.OriX;
    clear param.OriX;
end

if (~isfield(param,'InU'))
    s = median(abs(InX(:)));
    s = sqrt(s/r);
    if min(InX(:)) >= 0
        InU = rand(d,r)*s;
    else
        InU = rand(d,r)*s*2-s;
    end
    
else
    InU = param.InU;
end

if (~isfield(param,'InV'))
    InV = rand(r,n)*s;
else
    InV = param.InV;
end

if (~isfield(param,'k'))
    k = 3;
else
    k = param.k;
end


if (~isfield(param,'tol'))
    tol = 1.0e-5;
else
    tol = param.tol;
end

%%Initialize GMM parameters
tempX=InX(:);
R = initialization(tempX',k);
[~,label(1,:)] = max(R,[],2);
R = R(:,unique(label));
model.mu = zeros(1,k);
model.Sigma = rand(1,k);
model.alpha = 1.1*ones(1,r);
model.lambda = 10;
nk = sum(R,1);
model.weight = nk/size(R,1);
% llh = -inf(1,maxiter);
converged = false;
TempU = InU;
TempV = InV;


t = 1;
%%%%%%%%%%%%%%%%Initialized E Step %%%%%%%%%%%%%%%%%%%
[TempU, TempV] = expectationUV(InX,R,TempU,TempV,model);
TempX = TempU*TempV;
Error = InX(:)-TempX(:);
[R, llh(t)] = expectationR(Error',model,TempU,TempV);
%%%%%%%%%%%%%%%%Initialized E Step %%%%%%%%%%%%%%%%%%%

while ~converged && t < maxiter
    t = t+1;
    
    %%%%%%%%%%%%%%%% M Step %%%%%%%%%%%%%%%%%%%
    [model] = maximizationWS(Error',R,model);
    [model,k,R] = Ktuning(model,k,R);
    %%%%%%%%%%%%%%%% M Step %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%
    [TempU, TempV] = expectationUV(InX,R,TempU,TempV,model);
    TempX = TempU*TempV;
    Error = InX(:)-TempX(:);
    [R, llh(t)] = expectationR(Error',model,TempU,TempV);
    L1 = llh(t);
    %%%%%%%%%%%%%%%% E Step %%%%%%%%%%%%%%%%%%%
    
    if (~isfield(param,'k'))
        [~,label(:)] = max(R,[],2);
        u = unique(label);   % non-empty components
        if size(R,2) ~= size(u,2)
            R = R(:,u);   % remove empty components
        else
            converged = abs(llh(t)-llh(t-1)) < tol*abs(llh(t));
        end
        k = length(u);
    else
        converged = abs(llh(t)-llh(t-1)) < tol*abs(llh(t));
    end
    
end
OutU = TempU;
OutV = TempV;

disp(['The likelihood in this step is ',num2str(L1),';']);
disp(['There are ',num2str(k),' Gaussian noises mixed in data']);
disp(['with variances ',num2str(model.Sigma)]);
disp(['with weights ',num2str(model.weight),'.']);
disp(['Relative reconstruction error ', num2str(sum(sum(((OriX - TempU*TempV)).^2))/sum(sum((OriX).^2)))]);
disp(['L2 RMSE is ', num2str(sqrt(mean(Error.^2)))]);
disp(['L1 RMSE is ', num2str(mean(abs(Error)))]);

if converged
    fprintf('Converged in %d steps.\n',t-1);
else
    fprintf('Not converged in %d steps.\n',maxiter);
end