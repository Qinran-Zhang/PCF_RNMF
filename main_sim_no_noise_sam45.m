clc,clear,close all
%%
%Simulation data 
[ X, W, H, K, module_set]=create_simu_data(2);
sampleNum=size(W,1);
%I is the number of columns of each matrix in X
if iscell(X)
    I=zeros(1,length(X));
    for i=1:length(X)
        I(i)=size(X{i},2);
    end
else
    I=size(X,2);
end
%%
rand('seed',0)
exp_maxtimes=100;
%
module_by_NMF=cell(exp_maxtimes,1);
module_by_RNMF=cell(exp_maxtimes,1);
eucl_dist_by_NMF=zeros(exp_maxtimes,1);
eucl_dist_by_RNMF=zeros(exp_maxtimes,1);
caltime_by_NMF=zeros(exp_maxtimes,1);
caltime_by_RNMF=zeros(exp_maxtimes,1);
Ratio_by_NMF=zeros(exp_maxtimes,1);
Ratio_by_RNMF=zeros(exp_maxtimes,1);
maxiter=1000;
speak=0;
for exp_item=1:exp_maxtimes
    %NMF
    [W1, H1, eucl_dist_by_NMF(exp_item), ~, caltime_by_NMF(exp_item)] = NMF_re(X, K, 0,maxiter, speak);
    module_by_NMF{exp_item}=find_module(W1,H1,0);
    %RNMF
    [W2, H2, eucl_dist_by_RNMF(exp_item), ~, caltime_by_RNMF(exp_item)] = RNMF_re(X, K, 0,maxiter, 500, speak);
    module_by_RNMF{exp_item}=find_module(W2,H2,0);
    [ Order1,Ratio1, Ratio_by_NMF(exp_item) ]= eval_fun(module_by_NMF{exp_item},module_set);
    [ Order2,Ratio2, Ratio_by_RNMF(exp_item) ]= eval_fun(module_by_RNMF{exp_item},module_set);
end
disp(['Average error of NMF is ',num2str(mean(eucl_dist_by_NMF)),'Average error of RNMF is ',num2str(mean(eucl_dist_by_RNMF))]);
disp(['Average computational time of NMF is ',num2str(mean(caltime_by_NMF)),'Average computational time of NMF is ',num2str(mean(caltime_by_RNMF))]);
disp('Variance of error of NMF is ');
var(eucl_dist_by_NMF)
disp('Variance of error of RNMF is ');
var(eucl_dist_by_RNMF)
%%
%Similarity
disp(['Similarity of NMF is ',num2str(mean(Ratio_by_NMF)),'Similarity of RNMF is ',num2str(mean(Ratio_by_RNMF))]);
%%
%R value
[ P1,pmatrix1 ] = calPvalueSample( module_by_NMF, sampleNum);
[ P2,pmatrix2 ] = calPvalueSample( module_by_RNMF, sampleNum);
disp(['R value of NMF is ',num2str(P1),'R value of RNMF is ',num2str(P2)]);


% save 20180903_1859 X module_set module_by_NMF module_by_RNMF Ratio_by_NMF Ratio_by_RNMF eucl_dist_by_NMF eucl_dist_by_RNMF caltime_by_NMF caltime_by_RNMF 


