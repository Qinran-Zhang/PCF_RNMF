function [ Y ] = doublecol( X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[m,n]=size(X);
Y=zeros(m,2*n);
Y(:,1:2:end)=X.*(X>0);
Y(:,2:2:end)=-X.*(X<0);
end

