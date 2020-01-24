%不同搜索速度与经航速度之间的关系
clear all

%经航速度
Vse = 5/2;

%搜索概率
%数值积分
P = zeros(1,200);
k = zeros(1,200);
%离海面高度
h1 = -0;
%潜深(m)
h2 = 0;
%磁探仪作用距离
d = 1000;
W = 2*sqrt(d^2 - (h1+h2)^2);  
%距离目标
D0 = 30*10^3;%初始距离影响很大，应尽量减小
V0 = 4000;%假设的速度
Kd = 1;
T = 2*3600;

i = 1;
Max = 0;
n = [];
%k = [];
R2 = [];

for Vs = 25/2:1:280/2%搜索速度影响很大，应尽量提高
    t0 = D0/(Vse+Vs);
    R0 = Vse * t0;
    D = D0*Kd;%假设的距离
    t1 = D /(Vse + 10*20);
    R1 = Vse * t1;%假设的初始距离
    %k = [k,tan(asin(Vse/Vs))];
    k = tan(asin(Vse/Vs));
    n = [n,log(R0/R1)/(2*k*pi)];
    m = 0;
    while(1)
        m = m+1;
        R2 = R0*exp(k*2*m*pi);
        if R2 > Vse * (T-t0)  %超过搜索时间
            break
        end
    end

    if m>n(end)
        M = [m,i];
        %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R0),2*n(end)*pi,(2*m)*pi);       %#ok<DQUAD>
        %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R1),2*n(end)*pi,(2*n(end)+10)*pi);       %#ok<DQUAD>
        %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R0),2*0*pi,(2*0+16)*pi);       %#ok<DQUAD>
        P(i) = quad(@(x) ss_probability(x,k,Vse,W,R1),2*0*pi,(2*m+0)*pi);       %#ok<DQUAD>
        %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R1),2*n(end)*pi,(2*m)*pi);       %#ok<DQUAD>
        
        if Max < P(i)
            Max = P(i);
            V_m = Vs;
        end
    end
    i = i+1;
    if i<100
        R0 = D * Vse /(Vse + Vs);
        k = tan(asin(Vse/Vs));
        L = R0*exp(k*4.5*pi);
    end
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