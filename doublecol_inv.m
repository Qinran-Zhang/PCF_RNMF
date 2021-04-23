function [ X ] = doublecol_inv( Y ,para)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[m,n]=size(Y);
if nargin==1
    para=0;
end
if para==0
X=Y(:,1:2:end)-Y(:,2:2:end);
else
    Y1=Y(:,1:2:end);
    Y2=Y(:,2:2:end);
    flag=abs(Y1)>abs(Y2);
    Y1(~flag)=0;
    Y2(flag)=0;
    X=Y1-Y2;
end

