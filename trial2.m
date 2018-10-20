clear all; close all; clc;

%第二级优化，即云上已有东西时的优化。（State-transition Pricing）
%请于进行过第一级优化后使用
%请更改ALPHA_NEW, BETA_NEW, ROU_NEW, THETA_NEW, 以及t以调试

exp_times = 5;

c = 8;
d = 8;
s = 10;

x_C   = 4; %maximum number of components to be on the same cloud
x_D   = 4; %maximum number of data  sets to be on the same cloud

set_of_data = 3; %all-to-one, static, dynamic

Ti_min = 0; %10's power
Ti_max = 5;
Tinum  = 6;
Ttimes = 30;
Ti = ( logspace(Ti_min,Ti_max,Tinum) )';
Ti = floor(Ti);
T_set = Ttimes*Ti;
Til = length(Ti);
Tr = T_set(end);

% T = 1000;
% Ti_min = 10;
% Ti_inc = Ti_min;
% Ti_max = 50;
% Ti = (Ti_min:Ti_inc:Ti_max)';
% Til = length(Ti);

P = zeros(set_of_data,Til);

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
%     while (isnan(B(1,1)))
%         [ A,B ] = test( C,D,L,M,ALPHA,BETA,ROU,THETA,w,x_C,x_D );
%         disp('Error1 Occurred.')
%     end
    A2 = round(full(A)); B2 = round(full(B)); 

    for i = 1:Til
        tic;

        T = T_set(i);
        
        tT = (0:Ti(i):T)'; %set price changing moments
        if(tT(end) ~= T)
            tT = [tT; T]; %#ok<AGROW>
        end
        tTl = length(tT);
        
        j = 1;
        t = tT(j+1)-tT(j);

        p1 = price(A1,B1,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
        p2 = price(A2,B2,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
        A3 = A2; B3 = B2;
        p3 = price(A3,B3,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
        
        for j = 2:(tTl-1)
            t = tT(j+1)-tT(j);

            ALPHA = ALPHA_SET(:,tT(j));
            BETA  = BETA_SET(:,tT(j));
            ROU   = ROU_SET(:,tT(j));
            THETA = THETA_SET(:,tT(j));

            [ A,B ] = test2(A3,B3,t,C,D,L,M,ALPHA,BETA,ROU,THETA,w);
%             while (isnan(B(1,1)))
%                 [ A,B ] = test2(A3,B3,t,C,D,L,M,ALPHA,BETA,ROU,THETA,w);
%                 disp('Error2 Occurred.')
%             end
            
            r = trans_price( B3,round(full(B)),D,ROU,THETA );
            
            A3 = round(full(A)); B3 = round(full(B));
        
            p1 = p1 + price(A1,B1,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
            p2 = p2 + price(A2,B2,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t;
            p3 = p3 + price(A3,B3,C,D,L,M,ALPHA,BETA,ROU,THETA,w)*t + r;
        end

        P(:,i) = P(:,i) + [p1;p2;p3];

        toc;

    end

end

P = P/exp_times;

T_factor = Ti(end) ./ (Ti');
T_factor = repmat(T_factor,set_of_data,1);

P = P .* T_factor;

P = P / Ttimes / 10^Ti_max;

figure;
% plot(Ti,P(1,:),'v-k', ...
%      Ti,P(2,:),'s-b', ...
%      Ti,P(3,:),'^-r', ...
%      'LineWidth',2);
 
semilogx(Ti,P(1,:),'v-k', ...
         Ti,P(2,:),'s-b', ...
         Ti,P(3,:),'^-r', ...
         'LineWidth',2);

grid;
xlabel('Price Variety Time Interval t (hours)', ...
       'FontSize',20,'FontName','Calibri');
ylabel('Average Cloud Service Fee ($/hour)', ...
       'FontSize',20,'FontName','Calibri');

legend('Single-Cloud','Optimal-Static','Optimal-Cognitive', ...
       'Location','Best');

set(gca,'FontSize',15);

xlim([10^Ti_min,10^Ti_max]);



figure;
% plot(Ti,P(1,:),'v-k', ...
%      Ti,P(2,:),'s-b', ...
%      Ti,P(3,:),'^-r', ...
%      'LineWidth',2);
 
semilogx(                 ...
         Ti,P(3,:),'^-r', ...
         'LineWidth',2);

grid;
xlabel('Pricing Period (hours)', ...
       'FontSize',20,'FontName','Calibri');
ylabel('Total Price ($)', ...
       'FontSize',20,'FontName','Calibri');

legend('Optimal-Cognitive', ...
       'Location','Best');

set(gca,'FontSize',15);

xlim([10^Ti_min,10^Ti_max]);



