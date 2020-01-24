%不同搜索速度与经航速度之间的关系
clear all

%经航速度
Vse = 5*2;

%搜索概率
%数值积分
P = zeros(1,200);
k = zeros(1,200);
%离海面高度
h1 = -50;
%潜深(m)
h2 = 100;
%磁探仪作用距离
d = 1000;
W = 2*sqrt(d^2 - (h1+h2)^2); 
%散布
R = 1*10^3;
%s搜索时间
T= 12*3600;
%距离目标
D = 18 *10^3;
i = 1;
Max = 0;
for Vs = 5*2:1:20*2
    P(i) = 1-exp(-W*Vs*T/(pi*R*(R+Vse*T)));
    i = i+1;
    if Max < P(i-1)
        Max = P(i-1);
        V_m = Vs;
    end
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
figure(9)
plot(P(:))
V_m
Max