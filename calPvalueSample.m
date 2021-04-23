function [ P,pmatrix ] = calPvalueSample( module ,sampleNum)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
n=sampleNum;
K=size(module{1},1);
exp_times=length(module);
pmatrix=zeros(sampleNum,sampleNum);
for exp_item=1:exp_times
    cmat=false(sampleNum,sampleNum);
    for k=1:K
        module_one=module{exp_item}{k,1};
        for i=1:length(module_one)
            for j = 1:length(module_one)
                cmat(module_one(i),module_one(j))=true;
            end
        end
    end
    pmatrix=pmatrix+cmat;
end
pmatrix=pmatrix/exp_times;
p_nonzeros=nonzeros(pmatrix);
P=4/(n^2)*(sum((p_nonzeros-0.5).^2)+0.25*(n^2-length(p_nonzeros)));
end

