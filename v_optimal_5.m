%��ͬ�����ٶ��뾭���ٶ�֮��Ĺ�ϵ
clear all

%�����ٶ�
Vse = 8*0.5144444;
%��ʼɢ��
s0 = 1*10^3;
%��������
%��ֵ����
%P = zeros(1,40);
P = ones(1,180);
%�뺣��߶�
h1 = -0;
%Ǳ��(m)
h2 = 100;
%��̽�����þ���
d = 455;
W = 2*sqrt(d^2 - (h1+h2)^2);  
%W = d;
%����Ŀ��
D0 = 10*1.852*10^3;%��ʼ����Ӱ��ܴ�Ӧ������С
V0 = 0;%������ٶ�
Kd = 1;
T = 4*3600;  %����ʱ��

i = 1;
Max = 0;
n = [];
%k = [];

for Vs = 10:1:180 %�����ٶ�Ӱ��ܴ�Ӧ�������
    %R0 = Vse * t0;
    Vs1 = Vs*0.51444;%���㵽m
    D = D0*Kd;%����ľ���
    t0  = s0/(Vse);
    %t1 = D /(Vse + 300);
    t1 = (D-s0)/(Vse + Vs1);
    R0 = Vse * (t1+0); %����ĳ�ʼ����
    R1 = Vse * (t1+0); %����ĳ�ʼ����
%     %t1 = D /(Vse + 300);
    %k = [k,tan(asin(Vse/Vs))];
    k = tan(asin(Vse/Vs1));
    n = [n,log(R0/R1)/(2*k*pi)];
    m = 0;
    T1 = T-t1;
    M = floor((log(Vse * (T1)/R0)/(k*2*pi)))
    N = (log(Vse * (T1)/R0)/(k*2*pi))-M
    while(1)
        R2 = R0*exp(k*2*m*pi); %ת��mȦ֮�������ߵļ���
        for j = 1:M
            %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R0),2*n(end)*pi,(2*m)*pi);       %#ok<DQUAD>
            %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R1),2*n(end)*pi,(2*n(end)+10)*pi);       %#ok<DQUAD>
            %P(i) = P(i)*(1-quad(@(x) ss_probability(x,k,Vse,W,R0),2*j*pi,(2*j+2)*pi));       %#ok<DQUAD>
            %P(i) = P(i)+ss_probability_5(k,Vse,W,R1,t1+t0,0,N);    %�����ۼ�,P(i)��ʼ��Ϊ0
            P(i) = P(i)*(1-ss_probability_5(k,Vse,W,R1,t0,0,N));   
            %P(i) = quad(@(x) ss_probability(x,k,Vse,W,R1),2*n(end)*pi,(2*m)*pi);       %#ok<DQUAD>
            R1 = R1*exp(k*2*pi);
        end
        %P(i) = P(i)+ss_probability_5(k,Vse,W,R1,t1+t0,1,N);
        P(i) = P(i)*(1-ss_probability_5(k,Vse,W,R1,t0,1,N)); 
        R1 = R1*exp(k*2*N*pi);
        t2 = t1;
        t1 = R1/Vse;
        if t1 >= T1  %��������ʱ��
            break
        end
        %t0 = t1-t2
    end
    
    %P(i) = 1-(1-P(i))^4;
    P(i) = 1-P(i)^1;%4�����ķ��ָ���
    i = i+1;
    P
end
figure(2)
P4_180_far = P;
save('1_180_near.mat','P1_180_near');
plot(10:i+8,P(1:i-1))
%P
