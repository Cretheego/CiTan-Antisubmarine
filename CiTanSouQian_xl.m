%巡逻线长度(m)
W = 50*10^3; H = 50*10^3;
LD = 100*10^3;
%发现次数
%仿真次数
h = 1;
N = 1000/h;
find_cnt = zeros(1,N);

%搜潜兵力速度
V_s = 50*10^3;
W_m = 800;
%初始位置
S_l = [0,W_m/2];
%W_m是磁探仪的搜索宽度
%转弯的时间
t_zw = 5;
%周期T
T = 2* (LD/V_s)*3600 + 2*t_zw;
t1 = (LD/V_s)*3600;
t2 = t1 + t_zw;
t3 = 2* t1 + t_zw;
for n = 1:h:N
    %潜艇参数
    %初始位置，均匀分布
    %Q_l = [x1,y1];
    x_t = LD * rand(1);
    y_t = -5*10^3 - 2*10^3 * rand(1);
    %速度(节) 5-13kn
    V = (5+(13-5)*rand(1))*2;
   
    %航向垂直巡逻线
    H0 = pi/2;
    t = 0;
    while(1)
        %潜艇的位置
        x_t = x_t + V*cos(H0)*t;
        y_t = y_t + V*sin(H0)*t;
               
        %时刻t的磁探仪的位置
        t_T = t - floor(t/T)*T;
        x_st = V_s*t_T*(t_T>0 & t_T<=t1)/3600 + 10^8*(t_T>t1& t_T<=t2) + V_s*(t3-t_T)*(t_T>t2 & t_T<=t3)/3600 + 10^8*(t_T>t3& t_T<=T) ;
        y_st = (-W_m/2 )*(t_T>0 & t_T<=t2) + (W_m/2)*(t_T>t2 & t_T<=T);
        if (y_st > H)
            break;
        end
        %潜艇和探潜兵力之间的距离
        d_fq = sqrt((x_st - x_t)^2 + (y_st - y_t)^2);
        if d_fq <= W_m/2
            find_cnt(n) = find_cnt(n)+1;
            %Continue;
        end
        t = t+1;
        if y_t > 2*10^3
            break;
        end
    end
    n;
end
%发现概率
A = find_cnt>0;
find_pro = 1*sum(A(:))/N
