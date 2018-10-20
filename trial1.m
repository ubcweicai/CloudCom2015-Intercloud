clear all; close all; clc;

%第一级优化，即云上未有东西时的优化。（未考虑第一次上传价格）
%请更改C, D, ALPHA, BETA, ROU, THETA, L, M, 以及w以调试

exp_times = 20;

c = 8;
d = 8;
s = 10;

% x_C   = 8; %maximum number of components to be on the same cloud
x_C_set = ([1,2,4,6,8])'; %maximum number of components to be on the same cloud
% x_D   = 8; %maximum number of data  sets to be on the same cloud
x_D_set = (1:d)';

x_Cl = length(x_C_set);
x_Dl = length(x_D_set);

p2 = 0; 
p3 = zeros(x_Dl,x_Cl);

for li = 1:exp_times

    C     = randi([1,30],c,1);
    D     = randi([1024,10240],d,1);

    Lp = 0.1; L = rand([c,d]) < Lp;
    M = 0.1*rand(c);
    w = 0.00001;
    
    rand_dev = 0.1; 
    ALPHA = U(0.0081,0.04,rand_dev,s,1);%0.0319*rand(s,1)+0.0081;
    BETA  = U(0.000036,0.000056,rand_dev,s,1);%0.000020*rand(s,1)+0.000036;
    ROU   = U(0.0012,0.0021,rand_dev,s,1);%0.0009*rand(s,1)+0.0012;
    THETA = U(0.12,0.21,rand_dev,s,1);%0.09*rand(s,1)+0.12;
%     ALPHA = U(0.0081,0.04,s,1);%0.0319*rand(s,1)+0.0081;
%     BETA  = U(0.000036,0.000056,s,1);%0.000020*rand(s,1)+0.000036;
%     ROU   = U(0.0012,0.0021,s,1);%0.0009*rand(s,1)+0.0012;
%     THETA = U(0.12,0.21,s,1);%0.09*rand(s,1)+0.12;

    %all to one cloud
    p2_set = zeros(s,1);
    for k = 1:s
        A2 = zeros(c,s); A2(:,k) = 1;
        B2 = zeros(d,s); B2(:,k) = 1;
        p2_set(k) = price(A2,B2,C,D,L,M,ALPHA,BETA,ROU,THETA,w);
    end
    p2_min = find(p2_set == min(p2_set));
    A2 = zeros(c,s); A2(:,p2_min(1)) = 1;
    B2 = zeros(d,s); B2(:,p2_min(1)) = 1;
    p2 = p2 + p2_set(p2_min(1));

    for i = 1:x_Dl
        tic;
        x_D = x_D_set(i);
        
        for j = 1:x_Cl

            x_C = x_C_set(j);

%             %random allocation
%             A1 = rand_loc(c,s);
%             B1 = rand_loc(d,s);
%             p1(i) = p1(i) + price(A1,B1,C,D,L,M,ALPHA,BETA,ROU,THETA,w);

            %best case
            [ A,B ] = test( C,D,L,M,ALPHA,BETA,ROU,THETA,w,x_C,x_D );
            A3 = full(A);
            B3 = full(B);
            p3(i,j) = p3(i,j) + price(A3,B3,C,D,L,M,ALPHA,BETA,ROU,THETA,w);

        end
        
        toc;
    end
    
end

% p1 = p1/exp_times; 
p2 = p2/exp_times; 
p2 = repmat(p2,x_Dl,1);

p3 = p3/exp_times;

figure;
plot(x_D_set,p2     ,'.-k', ...
     x_D_set,p3(:,1),'^-r', ...
     x_D_set,p3(:,2),'o-g', ...
     x_D_set,p3(:,3),'*-b', ...
     x_D_set,p3(:,4),'s-c', ...
     x_D_set,p3(:,5),'v-m', ...
     'LineWidth',2 );

xlabel('Data Set Security Level (\Phi)', ...
       'FontSize',20,'FontName','Calibri');
ylabel('Cloud Service Fee ($)', ...
       'FontSize',20,'FontName','Calibri');
   
grid;
% legend('Random Allocation','All-to-One-Cloud','Optimal');
legend('Single-Cloud', ...
       '\Theta = 1','\Theta = 2','\Theta = 4','\Theta = 6','\Theta = 8', ...
       'Location','BestOutside');

set(gca,'FontSize',15);

xlim([x_D_set(1) x_D_set(end)]);




