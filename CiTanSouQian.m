clear all
%搜索区域(m)
W = 50*10^3; H = 50*10^3;
%发现次数
%仿真次数
h = 1;
N = 500/h;
find_cnt = zeros(1,N);

%搜潜兵力速度
V_s = 500*10^3;
W_m = 400;
%初始位置
S_l = [0,W_m/2];
%W_m是磁探仪的搜索宽度
%时刻t的位置
t_zw = 5;
%周期T
T = 2* W/V_s*3600 + 2*t_zw;
t1 = (W/V_s)*3600;
t2 = t1 + t_zw;
t3 = 2* t1 + t_zw;
num = 1;
for n = 1:h:N
    %潜艇参数
    %初始位置，均匀分布
    %Q_l = [x1,y1];
    Tot = 12*3600;
    x_t = W * rand(1);
    y_t = H * rand(1);
    x_st = zeros(num,Tot+1);
    y_st = zeros(num,Tot+1);
    %速度(节) 5-13kn
    V = (5+(13-5)*rand(1))*2;
    %航向，[0,2*pi]均匀分布
    H0 = 2*pi*rand(1);
    t = 0;
    while(1)
        x_t = x_t + V*cos(H0)*t;
        y_t = y_t + V*sin(H0)*t;
        %与边界发生碰撞
        if x_t < 0 || x_t > W
            x_t = x_t - V*cos(H0)*t;
            H0 = 2*pi - H0;
        end
        if y_t < 0 || y_t > H
            y_t = y_t - V*sin(H0)*t;
            H0 = pi - H0;
        end
        
        %时刻t的位置
        t_T = t - floor(t/T)*T;
        x_st(1,t+1) = V_s*t_T*(t_T>0 & t_T<=t1) + 0*(t_T>t1& t_T<=t2) + V_s*(t3-t_T)*(t_T>t2 & t_T<=t3) + 0*(t_T>t3& t_T<=T) ;
        y_st(1,t+1) = (W_m/2 + floor(t/T) *W_m)*(t_T>0 & t_T<=t2) + (3*W_m/2 + floor(t/T) *W_m)*(t_T>t2 & t_T<=T);
        if (y_st(1,t+1) > H)
            break;
        end
        %潜艇和探潜兵力之间的距离
        d_fq = sqrt((x_st(1,t+1) - x_t)^2 + (y_st(1,t+1) - y_t)^2);
        if d_fq <= W_m/2
            find_cnt(n) = find_cnt(n)+1;
            %Continue;
        end
        t = t+1;        
    end
    n;
end
%发现概率
figure(4)
plot(x_st(1,:),y_st(1,:))
A = find_cnt>0;
find_pro = 1*sum(A(:))/N
