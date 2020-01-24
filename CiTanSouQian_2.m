clear all
%搜索区域(m)
W = 50*10^3; H = 50*10^3;
%发现次数
%仿真次数
h = 1;
N = 500/h;
find_cnt = zeros(1,N);

%搜潜兵力速度
V_s = 50*10^3;
W_m = 800;
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
    x_t = zeros(num,Tot+2);
    y_t = zeros(num,Tot+2);
    x_t(1) = W * rand(1);
    y_t(1) = H * rand(1);
    x_st = zeros(num,Tot+1);
    y_st = zeros(num,Tot+1);
    %速度(节) 5-13kn
    V = (5+(13-5)*rand(1))*2;
    %航向，[0,2*pi]均匀分布
    H0 = 2*pi*rand(1);
    t = 0;
    while(1)
        x_t(t+2) = x_t(t+1) + V*cos(H0)*t;
        y_t(t+2) = y_t(t+1) + V*sin(H0)*t;               
        %与边界发生碰撞
        if x_t(t+2) < 0 || x_t(t+1) > W
            x_t(t+2) = x_t(t+2) - V*cos(H0)*t;
            H0 = 2*pi - H0;
        end
        if y_t(t+2) < 0 || y_t(t+2) > H
            y_t(t+2) = y_t(t+2) - V*sin(H0)*t;
            H0 = pi - H0;
        end
        
        %时刻t的位置
        t_T = t - floor(t/T)*T;
        x_st(1,t+1) = V_s*t_T*(t_T>0 & t_T<=t1) + 0*(t_T>t1& t_T<=t2) + V_s*(t3-t_T)*(t_T>t2 & t_T<=t3) + 0*(t_T>t3& t_T<=T) ;
        y_st(1,t+1) = (W_m/2 + floor(t/T) *W_m)*(t_T>0 & t_T<=t2) + (3*W_m/2 + floor(t/T) *W_m)*(t_T>t2 & t_T<=T);
        if (y_st(1,t+1) > H)
            break;
        end
        x_st(3,t+1) = (W_m/2 + floor(t/T) *W_m)*(t_T>0 & t_T<=t2) + (3*W_m/2 + floor(t/T) *W_m)*(t_T>t2 & t_T<=T); 
        y_st(3,t+1) = V_s*t_T*(t_T>0 & t_T<=t1) + 0*(t_T>t1& t_T<=t2) + V_s*(t3-t_T)*(t_T>t2 & t_T<=t3) + 0*(t_T>t3& t_T<=T) ;
        
        x_st(2,t+1) =  x_st(1,t+1);
        y_st(2,t+1) = (H-(W_m/2 + floor(t/T)) *W_m)*(t_T>0 & t_T<=t2) + (H-(3*W_m/2 + floor(t/T) *W_m))*(t_T>t2 & t_T<=T);
        
        x_st(4,t+1) = (W-(W_m/2 + floor(t/T) *W_m))*(t_T>0 & t_T<=t2) + (W-(3*W_m/2 + floor(t/T) *W_m))*(t_T>t2 & t_T<=T); 
        y_st(4,t+1) = y_st(3,t+1);
        
        %潜艇和探潜兵力之间的距离
        d_fq = sqrt((x_st(:,t+1) - repmat(x_t(t+1),num,1)).^2 + (y_st(:,t+1) - repmat(y_t(t+1),num,1)).^2);
        I = (d_fq <= W_m/2);
        if ~isempty(find(I(:)~=0,1))
            find_cnt(n) = find_cnt(n)+1;
            %Continue;
        end

%         d_fq = sqrt((x_st(1,t+1) - x_t(t+2))^2 + (y_st(1,t+1) - y_t(t+2))^2);
%         if d_fq <= W_m/2
%             find_cnt(n) = find_cnt(n)+1;
%             %Continue;
%         end
        t = t+1;
        if t > Tot 
            break;
        end 
    end
    n;
end
%发现概率
figure(2)
plot(x_t,y_t)
figure(4)
plot(x_st(3,:),y_st(3,:))
A = find_cnt>0;
find_pro = 1*sum(A(:))/N
