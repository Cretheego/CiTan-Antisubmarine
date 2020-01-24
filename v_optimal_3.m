%不同搜索速度与经航速度之间的关系
clear all

%经航速度
Vse = 5/2;
%初始散布
s0 = 1*10^3;
%搜索概率
%数值积分
P = ones(1,200);
k = zeros(1,200);
%离海面高度
h1 = -0;
%潜深(m)
h2 = 100;
%磁探仪作用距离
d = 1100;
W = 2*sqrt(d^2 - (h1+h2)^2);  
%距离目标
D0 = 10*10^3;%初始距离影响很大，应尽量减小
V0 = 4000;%假设的速度
Kd = 1;
T = 3*3600;  %搜索时间

i = 1;
Max = 0;
n = [];
%k = [];

for Vs = 20/2:1:340/2%搜索速度影响很大，应尽量提高
    t0 = D0/(Vse+1*Vs);
    R0 = Vse * t0;
    D = D0*Kd;%假设的距离
    t1 = D /(Vse + Vs);
    %t1 = D /(Vse + 300);
    R1 = Vse * t1; %假设的初始距离
    %k = [k,tan(asin(Vse/Vs))];
    k = tan(asin(Vse/Vs));
    n = [n,log(R0/R1)/(2*k*pi)];
    m = 1;
    while(1)
        m = m+1;
        R2 = R0*exp(k*2*m*pi);
        if R2 > Vse * (T-t1)  %超过搜索时间
            break
        end
    end

    if m>n(end)
        M = [m,i];
        n1 = n(end);
        for j = 0:m
            %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R0),2*n(end)*pi,(2*m)*pi);       %#ok<DQUAD>
            %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R1),2*n(end)*pi,(2*n(end)+10)*pi);       %#ok<DQUAD>
            %P(i) = P(i)*(1-quad(@(x) ss_probability(x,k,Vse,W,R0),2*j*pi,(2*j+2)*pi));       %#ok<DQUAD>
            %triplequad(@Fun,xm,xM,ym,yM,zm.zM)  %#ok<DTRIQD>
            P(i) = P(i)*(1-integral3(@(r,x,ss)ss_probability_3(r,x,ss,k,Vse,Vs,W,D,s0),0,inf,0,2*pi,0,2*pi));       %#ok<DTRIQD,DQUAD>
            %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R1),2*n(end)*pi,(2*m)*pi);       %#ok<DQUAD>
            R1 = R1*exp(k*2*pi);
        end
        P(i) = 1-P(i)^1;
%         if Max < P(i)
%             Max = P(i);
%             V_m = Vs;
%         end
    end
    i = i+1;
%     if i<100
%         R0 = D * Vse /(Vse + Vs);
%         k = tan(asin(Vse/Vs));
%         L = R0*exp(k*4.5*pi);
%     end
end
% fai =0.5*pi;
% u2 = 10;
% u1 = 2;
% Vse = (u2+u1)/2;
% W = 2.8;
% D =10;
% i=1;
% for Vs = 15:1:100*1
%     k(i) = tan(asin(Vse/Vs));
%     R0 = D * Vse /(Vse + Vs);
%     u11 = Vse - W*Vse/(2*R0);
%     if u11 < u1
%         u11 = u1;
%     end
%            
%     u22 = Vse + W*Vse/(2*R0);
%     if u22>u2
%         u22 = u2;
%     end
%     
%     P(i) = W*Vse*(1-exp(-k(i)*fai))/(2*k(i)*pi*(u22-u11)*R0);
%     i = i+1;
%     if Max < P(i-1)
%         Max = P(i-1);
%         V_m = Vs;
%     end
%     if i<100
%         R0 = D * Vse /(Vse + Vs);
%         L = R0*exp(k(i)*4.5*pi);
%     end
% end
figure(8)
plot(P(:))
figure(9)
subplot(211)
plot(n)
subplot(212)
plot(R2)
%V_m
Max
R2
k
i