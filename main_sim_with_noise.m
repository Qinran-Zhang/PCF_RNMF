clc,clear,close all
%%
%Simulation data 
[ X, W, H, K, module_set ]=create_simu_data(7);
%I is the number of columns of each matrix in X
X=add_noise(X,3,0.1);
if iscell(X)
    I=zeros(1,length(X));
    for i=1:length(X)
        I(i)=size(X{i},2);
    end
else
    I=size(X,2);
end
load('sim_data4');

%%
rand('seed',0)
exp_maxtimes=100;
lambda_all=[0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1];
module_by_NMF=cell(exp_maxtimes,length(lambda_all));
module_by_RNMF=module_by_NMF;
eucl_dist_by_NMF=zeros(exp_maxtimes,length(lambda_all));
eucl_dist_by_RNMF=eucl_dist_by_NMF;
caltime_by_NMF=eucl_dist_by_NMF;
caltime_by_RNMF=eucl_dist_by_NMF;
Ratio_by_NMF=eucl_dist_by_NMF;
Ratio_by_RNMF=eucl_dist_by_NMF;
BIC_by_NMF=eucl_dist_by_NMF;
BIC_by_RNMF=eucl_dist_by_NMF;
tol=1e-10;
maxiter=2000;
speak=1;
% lambda_all=[0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1];
%lambda_all=0.001;
for i=1:length(lambda_all)
    i
    lambda=lambda_all(i);
    for exp_item=1:exp_maxtimes
        %NMF
        [W1, H1, eucl_dist_by_NMF(exp_item,i), ~, caltime_by_NMF(exp_item,i)] = NMF_re(X, K, lambda, maxiter, speak);
        module_by_NMF{exp_item,i}=find_module(W1,H1,0);
        %RNMF
        [W2, H2, eucl_dist_by_RNMF(exp_item,i), ~, caltime_by_RNMF(exp_item,i)] = RNMF_re(X, K, lambda, maxiter, 500, speak);
        module_by_RNMF{exp_item,i}=find_module(W2,H2,0);
        sum(sum(W1>tol))
        sum(sum(W2>tol))
        [ Order1,Ratio1, Ratio_by_NMF(exp_item,i) ]= eval_fun(module_by_NMF{exp_item,i},module_set);
        [ Order2,Ratio2, Ratio_by_RNMF(exp_item,i) ]= eval_fun(module_by_RNMF{exp_item,i},module_set);
        BIC_by_NMF(exp_item,i)=2*log(eucl_dist_by_NMF(exp_item,i))+2*log(sum(sum(W1>tol)));
        BIC_by_RNMF(exp_item,i)=2*log(eucl_dist_by_RNMF(exp_item,i))+2*log(sum(sum(W2>tol)));
    end
end
disp('Average error of NMF is ');
mean(mean(eucl_dist_by_NMF))
disp('Average error of RNMF is ');
mean(mean(eucl_dist_by_RNMF))
disp('Variance of error of NMF is ');
var(reshape(eucl_dist_by_NMF,nnz(eucl_dist_by_NMF),1))
disp('Variance of error of RNMF is ');
var(reshape(eucl_dist_by_RNMF,nnz(eucl_dist_by_RNMF),1))

disp('Average computational time of NMF is ');
mean(mean(caltime_by_NMF))
disp('Average computational time of RNMF is ');
mean(mean(caltime_by_RNMF))
%%
disp('Similarity of NMF is ');
mean(mean(Ratio_by_NMF))
disp('Similarity of RNMF is ');
mean(mean(Ratio_by_RNMF))
%%
%R value
P1=zeros(1,length(lambda_all));
P2=P1;
for i=1:length(P1)
    P1(i)=calPvalue( module_by_NMF(:,i), I);
    P2(i)=calPvalue( module_by_RNMF(:,i), I);
end
disp('R value of NMF is ')
mean(P1)
disp('R value of RNMF is ')
mean(P2)

% save sim_data4


