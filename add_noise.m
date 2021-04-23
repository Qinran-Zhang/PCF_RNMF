function [ X ] = add_noise( X, para, sig )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin<2
    para=1;
end
if para==1
    for i=1:length(X)
        X{i}=X{i}+sig*randn(size(X{i}));
        X{i}(X{i}<0)=0;
    end
elseif para==2
    X{1}(:,1)=X{1}(:,1)+sig*randn(size(X{1},1),1);
    X{1}(X{1}<0)=0;
elseif para==3
    for i=1:length(X)
        m=size(X{i},2);
        index=randperm(m);
        index=index(1:round(m/50));
        X{i}(:,index)=X{i}(:,index)+sig*randn(size(X{i},1),round(m/50));
        X{i}(X{i}<0)=0;
    end
end

end

