function [W, H, eucl_dist, iter, caltime] =RNMF_re(X, K, lambda, maxiter, print_iter, speak)
%2018.11.17最新修改，第一版，完全按照论文中方法
narginchk(2,Inf);
if nargin<6
    speak=true;
end
if nargin<4
    maxiter=1000;
    print_iter=maxiter/5;
end
if nargin<3
    lambda=0;
end
% print_iter = 100;

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

R=rank(X{1})
Omega=rand(R,m);
Y=cell(1,I);
Q=cell(1,I);
XR=cell(1,I);
HR=cell(1,I);
for i=1:I
    Y{i}=Omega*X{i};
    [Q{i},~]=qr(Y{i}',0);
    XR{i}=doublecol(X{i}*Q{i});
    HR{i}=rand(K,size(XR{i},2));
%     HR{i}=doublecol(H{i}*Q{i});
end
% Y=Omega*X;
% [Q,~]=qr(Y',0);
% XX=doublecol(X*Q);
% HH=doublecol(H*Q);
% HH = HH.*(W'*XX)./((W'*W)*HH+eps);
% W = W.*(HH*XX')'./(W*(HH*HH')+eps);
% H = doublecol_inv(HH)*Q';
% H(H<0) = 0;


% use W*H to test for convergence
% Xr_old=W*H;

t1=clock;
for iter=1:maxiter
    % 更新各个H
    for i=1:I
        HR{i} = HR{i}.*(W'*XR{i})./((W'*W)*HR{i}+eps);
%         H{i} = doublecol_inv(HR{i})*Q{i}';
%         HR{i} = doublecol(H{i}*Q{i});
    end
    % 更新W
    HXt = zeros(K,m);
    HHt = zeros(K,K);
    for i=1:I
        HXt = HXt+HR{i}*XR{i}';
        HHt = HHt+HR{i}*HR{i}';
    end
    W = W.*(HXt)'./(W*(HHt)+lambda*W+eps);


    if (rem(iter,print_iter)==0)
        for i=1:I
            H{i} = doublecol_inv(HR{i})*Q{i}';
            H{i}(H{i}<0)=0;
            HR{i} = doublecol(H{i}*Q{i});
        end
%         Xr = W*H;
%         diff = sum(sum(abs(Xr_old-Xr)));
%         Xr_old = Xr;
        if speak
            eucl_dist = nmf_euclidean_dist(X,W,H,lambda);
            %         errorx= mean(mean(abs(X-W*H)))/mean(mean(X));
            disp(['Iter = ',int2str(iter),...
                ', eucl dist ' num2str(eucl_dist)])
        end
    end
end
t2=clock;
caltime = etime(t2, t1);
% eucl_dist = nmf_euclidean_dist(X,W,H,lambda)
Wm = max(W);
W = W*diag(1./Wm);
for i=1:I
    H{i}=diag(Wm)*H{i};
end
eucl_dist = nmf_euclidean_dist(X,W,H,lambda);




function err = nmf_euclidean_dist(X,W,H,lambda)
I=length(X);
err = 0;
for i=1:I
    err = err+norm(X{i}-W*H{i},'fro')^2;
end
err = err+lambda*norm(W,'fro')^2;


