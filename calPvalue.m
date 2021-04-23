function [ P,pmatrix ] = calPvalue( module ,I)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
J=[0 cumsum(I)];
n=sum(I);
K=size(module{1},1);
exp_times=length(module);
pmatrix=zeros(n,n);
for exp_item=1:exp_times
    cmat=false(n,n);
    for kk=1:K
        module_one=module{exp_item}(kk,2:end);
        for i=1:length(I)
            for j = 1:length(I)
                %module{exp_item,}
                for k=1:length(module_one{i})
                    for l=1:length(module_one{j})
                        cmat(J(i)+module_one{i}(k),J(j)+module_one{j}(l))=true;
                    end
                end
            end
        end
    end
    pmatrix=pmatrix+cmat;
end
pmatrix=pmatrix/exp_times;
p_nonzeros=nonzeros(pmatrix);
P=4/(n^2)*(sum((p_nonzeros-0.5).^2)+0.25*(n^2-length(p_nonzeros)));
end

