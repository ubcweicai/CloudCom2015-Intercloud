clear all; clc; %close all;

%第二级优化，即云上已有东西时的优化。（State-transition Pricing）
%请于进行过第一级优化后使用
%请更改ALPHA_NEW, BETA_NEW, ROU_NEW, THETA_NEW, 以及t以调试

exp_times = 25;

c = 8;
d = 8;
s = 10;

x_C   = 4; %maximum number of components to be on the same cloud
x_D   = 4; %maximum number of data  sets to be on the same cloud

T  = 50;
Ti = 1;
Tr = T/Ti;
set_of_data = 3; %all-to-one, static, dynamic

tT = (0:Ti:T)';
if(tT(end) ~= T)
    tT = [tT; T]; 
end
tTl = length(tT);

P = zeros(set_of_data,tTl);

for li = 1:exp_times
    tic;

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

    ALPHA_SET = U(0.0081,0.04,rand_dev,s,Tr);%0.0319*rand(s,T)+0.0081;
    BETA_SET  = U(0.000036,0.000056,rand_dev,s,Tr);%0.000020*rand(s,T)+0.000036;
    ROU_SET   = U(0.0012,0.0021,rand_dev,s,Tr);%0.0009*rand(s,T)+0.0012;
    THETA_SET = U(0.12,0.21,rand_dev,s,Tr);%0.09*rand(s,T)+0.12;

    %all to one cloud
    AB1_set = zeros(s,1);
    for k = 1:s
        A = zeros(c,s); A(:,k) = 1;
        B = zeros(d,s); B(:,k) = 1;
        AB1_set(k) = price(A,B,C,D,L,M,ALPHA,BETA,ROU,THETA,w);
    end
    p2_min = find(AB1_set == min(AB1_set));
    A1 = zeros(c,s); A1(:,p2_min(1)) = 1;
    B1 = zeros(d,s); B1(:,p2_min(1)) = 1;

    %best case
    [ A,B ] = test( C,D,L,M,ALPHA,BETA,ROU,THETA,w,x_C,x_D );
    A2 = round(full(A));
    B2 = round(full(B));
    
    j = 1;
    t = tT(j+1)-tT(j);
    
    p1 = price(A1,B1,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
    p2 = price(A2,B2,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
    A3 = A2;  B3 = B2;
    p3 = price(A3,B3,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;

    P(:,j+1) = P(:,j+1) + [p1;p2;p3];
    
    p1 = 0; p2 = 0; p3 = 0;
    for j = 2:(tTl-1)

        t = tT(j+1)-tT(j);

        ALPHA = ALPHA_SET(:,floor(tT(j)/Ti));
        BETA  = BETA_SET(:,floor(tT(j)/Ti));
        ROU   = ROU_SET(:,floor(tT(j)/Ti));
        THETA = THETA_SET(:,floor(tT(j)/Ti));

        [ A,B ] = test2(A3,B3,t,C,D,L,M,ALPHA,BETA,ROU,THETA,w);
        
        r = trans_price( B3,round(full(B)),D,ROU,THETA );
        
        A3 = round(full(A)); B3 = round(full(B));
        
        p1 = price(A1,B1,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
        p2 = price(A2,B2,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
        p3 = price(A3,B3,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t + r;

%         p1 = p1 + price(A1,B1,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
%         p2 = p2 + price(A2,B2,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
%         p3 = p3 + price(A3,B3,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;

        P(:,j+1) = P(:,j+1) + [p1;p2;p3];
    
    end
    
    toc;
end

P = P/exp_times;
    
figure;
plot(tT,P(1,:),'v-k', ...
     tT,P(2,:),'s-b', ...
     tT,P(3,:),'^-r', ...
     'LineWidth',2);

grid;
xlabel('Time Elapsed (hours)','FontSize',20,'FontName','Calibri');
ylabel('Instant Cloud Service Fee ($)','FontSize',20,'FontName','Calibri');

legend('Single-Cloud','Optimal-Static','Optimal-Cognitive', ...
       'Location','Best');

set(gca,'FontSize',15);




