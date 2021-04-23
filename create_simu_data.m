function [ X, W, H, K, module_set ] =create_simu_data( para)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if para==1
    W=[1 0;1 0;0 1];
    H=[1 1 0 0;0 0 1 1];
    X=W*H;
    K=2;
elseif para==2
    w=ones(1,10);
    h1=ones(1,30);  h2=ones(1,40); h3=ones(1,50);
    
    W=zeros(45,4);    % 5
    H1=zeros(4,130);  % 10
    H2=zeros(4,170);  % 10
    H3=zeros(4,215);  % 15
    K=4;
    module_set=cell(K,4);
    W(1:10,1)=w';  W(11:20,2)=w'; W(21:30,3)=w'; W(31:40,4)=w';
    module_set{1,1}=1:10; module_set{2,1}=11:20; module_set{3,1}=21:30; module_set{4,1}=31:40;
    H1(1,1:30)=h1; H1(2,31:60)=h1; H1(3,61:90)=h1; H1(4,91:120)=h1;
    module_set{1,2}=1:30; module_set{2,2}=31:60; module_set{3,2}=61:90; module_set{4,2}=91:120;
    H2(1,1:40)=h2; H2(2,41:80)=h2; H2(3,81:120)=h2; % H2(4,121:160)=h2;
    module_set{1,3}=1:40; module_set{2,3}=41:80; module_set{3,3}=81:120; %module_set{4,3}=121:160;
    H3(1,1:50)=h3; H3(2,51:100)=h3; % H3(3,101:150)=h3;
    H3(4,151:200)=h3;
    module_set{1,4}=1:50; module_set{2,4}=51:100; %module_set{3,4}=81:120; 
    module_set{4,4}=151:200;

    % Original matrix
    X1=W*H1; X2=W*H2; X3=W*H3;
    X={X1,X2,X3};
    H={H1,H2,H3};

    
elseif para==3
    W=rand(100,10);
    H=rand(10,1000);
    X=W*H;
    K=5;
elseif para==4
    W=zeros(100,10);
    H=zeros(10,100);
    for i=1:10
        W((i-1)*10+1:i*10,i)=ones(10,1);
        H(i,(i-1)*10+1:i*10)=ones(1,10);
    end
    X=W*H;
    K=10;
elseif para==5
    n=100;
    W=zeros(n,10);
    H=zeros(10,n);
    for i=1:10
        W((i-1)*n+1:i*n,i)=ones(n,1);
        H(i,(i-1)*n+1:i*n)=ones(1,n);
    end
    X=W*H;
    K=10;
elseif para==6
    W=zeros(100,10);
    H=zeros(10,1000);
    for i=1:10
        W((i-1)*10+1:i*10,i)=ones(10,1);
        H(i,(i-1)*100+1:i*100)=ones(1,100);
    end
    X=W*H;
    K=10;
elseif para==7
    w=ones(1,10);
    h=ones(1,100);
    K=10;
    I=110;
    J=[1100,1200,1300];
    W=zeros(110,10);
    H{1}=zeros(10,1100);
    H{2}=zeros(10,1200);
    H{3}=zeros(10,1300);
    module_set=cell(K,4);
    for i=1:10
        W(10*i-9:10*i,i)=w;
        module_set{i,1}=10*i-9:10*i;
        H{1}(i,100*i-99:100*i)=h;
        module_set{i,2}=100*i-99:100*i;
        module_set{i,3}=[];
        module_set{i,4}=[];
    end
    for i=1:2:5
        H{2}(i,100*i-99:100*i)=h;
        module_set{i,3}=100*i-99:100*i;
        H{3}(i,100*i-99:100*i)=h;
        module_set{i,4}=100*i-99:100*i;
    end
    for i=7:2:10
        H{2}(i,100*i-99:100*i)=h;
        module_set{i,3}=100*i-99:100*i;
    end
    for i=8:2:10
        H{3}(i,100*i-99:100*i)=h;
        module_set{i,4}=100*i-99:100*i;
    end
    X=cell(1,3);
    for i=1:3
        X{i}=W*H{i};
    end
elseif para==8
    W=[1 0;1 0;0 1;0 1];
    H=W';
    X=W*H;
    K=2;
end
end

