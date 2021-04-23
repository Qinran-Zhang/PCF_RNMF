% function [ Co_module, Subpattern1, Subpattern2, Subpattern3]=TriNMF_module2(X1, X2, X3, W, H1, H2, H3, tt) 
%
function module = find_module(W,H,t)
% % Just cover the Subpattern1 and Subpattern2 output. 
%
% Compute the mean(meadia) and std in columns of W and rows in H1, H2 to determine
% the module member and output the Co-module based on W and H1, H2 matrices.
% Co_module is the list of sample and SNPs, Genes in a Co-module.
% Subpattern1, Subpattern2 is the corresponding sub-matrices for each Co-module.
% 
% Output rule 2
n = size(W,1);
K = size(W,2);
min_threshold=0.2;
if ~iscell(H)
    module=cell(K,2);
    for i=1:K
        %模块中包含的样本
        threshold=graythresh(W(:,i));
        module{i,1}=find(W(:,i)>max(threshold,min_threshold))';
        %模块中包含的基因
        threshold=graythresh(H(i,:));
        module{i,2}=find(H(i,:)>max(threshold,min_threshold));
    end
else
    I=length(H);
    module=cell(K,I+1);
    for i=1:K
        threshold=graythresh(W(:,i));
        module{i,1}=find(W(:,i)>max(threshold,min_threshold))';
        for j=1:I
            threshold=graythresh(H{j}(i,:));
            module{i,j+1}=find(H{j}(i,:)>max(threshold,min_threshold));
        end
    end
    
end






% % MW =mean(W,1);     MH =mean(H,2);
% MW =median(W,1);   MH1 =median(H1,2); MH2 =median(H2,2); MH3 =median(H3,2);
% 
% % VW =std(W,0,1);    VH1 =std(H1,0,2);  VH2 =std(H2,0,2); VH3 =std(H3,0,2);
% VW =mad(W,1,1);    VH1 =mad(H1,1,2);  VH2 =mad(H2,1,2); VH3 =mad(H3,1,2);  %mad: Mean or median absolute deviation
% 
% % Co-Module
% for i=1:K
%    c1=find(H1(i,:)> MH1(i) + tt1*VH1(i));
%    module1{i,1}=c1'; 
% 
%    c2=find(H2(i,:)> MH2(i) + tt2*VH2(i));
%    module2{i,1}=c2'; 
%    
%    c3=find(H3(i,:)> MH3(i) + tt3*VH3(i));
%    module3{i,1}=c3'; 
%    
%    r=find(W(:,i)> MW(i) + tt0*VW(i));
%    Co_module{i,1}=r'; Co_module{i,2}=c1'; Co_module{i,3}=c2'; Co_module{i,4}=c3';

%    Subpattern1{i}=X1(r,c1');
%    Subpattern2{i}=X2(r,c2');
%    Subpattern3{i}=X2(r,c3');
end


