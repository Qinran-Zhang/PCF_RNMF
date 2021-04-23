function [W, H, D, iter, caltime] = NMF_re(X, K, lambda, maxiter, speak)
narginchk(2,Inf);
if nargin<5
    speak=true;
end
if nargin<4
    maxiter=1000;
end
if nargin<3
    lambda=0;
end
print_iter = 500;

if ~iscell(X)
    X={X};
end

I=length(X);
[m,n]=size(X{1});
for i=1:I
    n(i)=size(X{i},2);
end
W=rand(m,K);
H=cell(1,I);
for i=1:I
    H{i}=rand(K,n(i));
end


t1=clock;
for iter=1:maxiter
    % Euclidean multiplicative method
    for i=1:I
        H{i} = H{i}.*(W'*X{i})./((W'*W)*H{i}+eps);
    end
    HXt = H{1}*X{1}';
    HHt = H{1}*H{1}';
    for i=2:I
        HXt = HXt+H{i}*X{i}';
        HHt = HHt+H{i}*H{i}';
    end
    W = W.*(HXt)'./(W*(HHt)+lambda*W+eps);
    if (rem(iter,print_iter)==0) && speak,
        eucl_dist = nmf_euclidean_dist(X,W,H,lambda);
        errorxz=0;
        errorxm=0;
        for i=1:I
            errorxz=errorxz+mean(mean(abs(X{i}-W*H{i})));
            errorxm=errorxm+mean(mean(X{i}));
        end
        errorx = errorxz/errorxm;
        disp(['Iter = ',int2str(iter),...
            ', relative error = ',num2str(errorx),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
end
t2=clock;
caltime = etime(t2, t1);
% D = nmf_euclidean_dist(X,W,H,lambda)
Wm = max(W);
W = W*diag(1./Wm);
for i=1:I
    H{i} = diag(Wm)*H{i};
end
D = nmf_euclidean_dist(X,W,H,lambda);




function err = nmf_euclidean_dist(X,W,H,lambda)
I=length(X);
err = 0;
for i=1:I
    err = err+norm(X{i}-W*H{i},'fro')^2;
%     err = err+sum(sum((X{i}-W*H{i}).^2));
end
err = err+lambda*norm(W,'fro')^2;











