function [ Order,Ratio, Ratio_totle ] = eval_fun( Co_module,Comodule )
%计算两个模块的相似性
%Order为两者中模块的对应关系
K=size(Co_module,1);
Order=ones(K,1);
Ratio=ones(K,4);
numerator=0;
denominator=0;
for i=1:K
    num_union=length(intersect(Co_module{i,1},Comodule{1,1}));
    for j=2:K
        num_union1=length(intersect(Co_module{i,1},Comodule{j,1}));
        if num_union1>num_union
            num_union=num_union1;
            Order(i)=j;
        end
    end
    Ratio(i,1)=num_union/length(union(Co_module{i,1},Comodule{Order(i),1}));
    numerator=numerator+num_union;
    denominator=denominator+length(union(Co_module{i,1},Comodule{Order(i),1}));
    for j=2:4
        numerator=numerator+length(intersect(Co_module{i,j},Comodule{Order(i),j}));
        denominator=denominator+length(union(Co_module{i,j},Comodule{Order(i),j}));
        if(isempty(union(Co_module{i,j},Comodule{Order(i),j})))
            Ratio(i,j)=1;
        else
            Ratio(i,j)=length(intersect(Co_module{i,j},Comodule{Order(i),j}))/length(union(Co_module{i,j},Comodule{Order(i),j}));
        end
    end
end
Ratio_totle=numerator/denominator;
end

