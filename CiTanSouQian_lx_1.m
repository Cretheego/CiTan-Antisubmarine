%螺旋搜索
%初始点(m)
x0 = 0; y0 = 0;
%搜索者距离
D = 10*10^3;
%发现次数
%仿真次数
h = 1;
N = 100/h;

%搜潜兵力速度
Vs = 10;
%经航速度
Vsem = 5/2;
%螺旋频率
k = tan(asin(Vsem/Vs));
%W_m是磁探仪的搜索宽度
W_m = 1500;
T = 6*3600;
num = 4;
Tot = 0;
N1 = 100;
P = ones(N,N1); 
for n = 1:N
    %潜艇参数
    %初始位置，均匀分布
    %在散布圆上均匀分布
    %速度(节) 5-13kn
    %Vse = (3+2*rand(1))/2;
    %航向，[0,2*pi]均匀分布
    H0 = 2*pi*n/N;
    N1 =100;    
    
    for m=1:N1
        x_t = zeros(1,T+1);
        y_t = zeros(1,T+1);
        x_st = zeros(num,T+1);
        y_st = zeros(num,T+1);
        B = Vsem*sqrt(2/pi);
        Vse = raylrnd(B,1,1);
        %螺旋频率
        k = tan(asin(Vsem/Vs));
        t = 1;
        t1 = D /(Vse + Vs);
        R1 = Vse * t1; 
        %R0 = 0;
        alpha = (2*pi/num:2*pi/num:2*pi)';
        j = 1;
        find_cnt = zeros(1,N1);
        T0 = 0;
        while(1)%每次循环
            x_t(t) = x_t(1) + Vse*cos(H0)*t;
            y_t(t) = y_t(1) + Vse*sin(H0)*t;

            %时刻t磁探仪的位置
            %find_cnt = zeros(1,N1);
            if t<t1
                x_st(:,t) = (D-Vs*(t)).*cos(alpha(:,1));
                y_st(:,t) = (D-Vs*(t)).*sin(alpha(:,1));
                %潜艇和探潜兵力之间的距离
                d_fq = sqrt((x_st(:,t) - repmat(x_t(t),num,1)).^2 + (y_st(:,t) - repmat(y_t(t),num,1)).^2);
                I = (d_fq <= W_m/2);
                if ~isempty(find(I(:)~=0,1))
                    find_cnt(j) = find_cnt(j)+1;
                    %Continue;
                end
            else
                if t ==t1
                    P(n,m) = P(n,m)*(1-find_cnt(1)/t1);
                end
                Rou = Vse * t;
                T0 = T0+1;
                sita = log(Rou/R1)/k+alpha;
                if log(Rou/R1)/k>2*j*pi
                    P(n,m) = P(n,m)*(1-find_cnt(j+1)/T0);
                    j = j+1;
                   % t1 = t;
                    T0 = 0;
                end
                x_st(:,t) = Rou.*cos(sita(:,1));
                y_st(:,t) = Rou.*sin(sita(:,1));
                %潜艇和探潜兵力之间的距离
                d_fq = sqrt((x_st(:,t) - repmat(x_t(t),num,1)).^2 + (y_st(:,t) - repmat(y_t(t),num,1)).^2);
                I = (d_fq <= W_m/2);
                if ~isempty(find(I(:)~=0,1))
                    find_cnt(j+1) = find_cnt(j+1)+1;
                    %T_end = t;
                    %Continue;
                end
            end
            t = t+1;
            if t > T 
                break;
            end 
        end
        n
        %Tot = Tot+T_end;
    end
end
%发现概率
A = find_cnt;
figure(1)
plot(x_t,y_t,x_st,y_st)
1-mean(mean(P))
