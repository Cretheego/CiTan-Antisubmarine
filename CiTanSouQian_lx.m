%螺旋搜索
%初始点(m)
x0 = 0; y0 = 0;
%散布距离
L = 2*10^3;
%发现次数
%仿真次数
h = 1;
N = 500/h;
find_cnt = zeros(1,N);
%搜潜兵力速度
V_s = 500*10^3;
%螺旋频率
A0 = L;
w = 1;
%W_m是磁探仪的搜索宽度
W_m = 400;
%搜索时间
T = 12*3600;
%外圈大小
V_max = 14;
R_max = T*V_max*2;
%搜索兵力数量
num = 1;
for n = 1:h:N
    %潜艇参数
    %初始位置，均匀分布
    %在散布圆上均匀分布
    H = 2*pi*rand(1);
    
    x_t = zeros(1,T+1);
    y_t = zeros(1,T+1);
    x_st = zeros(num,T+1);
    y_st = zeros(num,T+1);
    x_t(1) = L * cos(H);
    y_t(1) = L * sin(H);
    x_st(:,1) = R_max*cos(2*pi*(1/num:1/num:1))';
    y_st(:,1) = R_max*sin(2*pi*(1/num:1/num:1))';   
    A0 = R_max;
    %速度(节) 5-13kn
    V = (5+(13-5)*rand(1))*2;
    %航向，[0,2*pi]均匀分布
    H0 = 2*pi*rand(1);
    t = 0;
    sita = 0;
    %减小的间距
    A1 = (R_max - L) / T;
    while(1)
        x_t(t+1) = x_t(t+1) + V*cos(H0)*t;
        y_t(t+1) = y_t(t+1) + V*sin(H0)*t;
               
        %时刻磁探仪t的位置
        w = V_s / A0;

        A0 = A0 - 1*A1;
        if (mod(t,10)==0)
            x_st(:,t+1) =  A0 *cos(w*t + sita);
            y_st(:,t+1) =  A0 *sin(w*t + sita);
        end
        %潜艇和探潜兵力之间的距离
        d_fq = sqrt((x_st(t+1) - x_t(t+1))^2 + (y_st(t+1) - y_t(t+1))^2);
        if d_fq <= W_m/2
            find_cnt(n) = find_cnt(n)+1;
            %Continue;
        end
        t = t+1; 
        if t > T 
            break;
        end 
    end
    
    n;
end
A0
%发现概率
A = find_cnt>0;
figure(1)
plot(x_t,y_t,x_st,y_st)
find_pro = 1*sum(A(:))/N
