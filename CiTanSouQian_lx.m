%��������
%��ʼ��(m)
x0 = 0; y0 = 0;
%ɢ������
L = 2*10^3;
%���ִ���
%�������
h = 1;
N = 500/h;
find_cnt = zeros(1,N);
%��Ǳ�����ٶ�
V_s = 500*10^3;
%����Ƶ��
A0 = L;
w = 1;
%W_m�Ǵ�̽�ǵ��������
W_m = 400;
%����ʱ��
T = 12*3600;
%��Ȧ��С
V_max = 14;
R_max = T*V_max*2;
%������������
num = 1;
for n = 1:h:N
    %Ǳͧ����
    %��ʼλ�ã����ȷֲ�
    %��ɢ��Բ�Ͼ��ȷֲ�
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
    %�ٶ�(��) 5-13kn
    V = (5+(13-5)*rand(1))*2;
    %����[0,2*pi]���ȷֲ�
    H0 = 2*pi*rand(1);
    t = 0;
    sita = 0;
    %��С�ļ��
    A1 = (R_max - L) / T;
    while(1)
        x_t(t+1) = x_t(t+1) + V*cos(H0)*t;
        y_t(t+1) = y_t(t+1) + V*sin(H0)*t;
               
        %ʱ�̴�̽��t��λ��
        w = V_s / A0;

        A0 = A0 - 1*A1;
        if (mod(t,10)==0)
            x_st(:,t+1) =  A0 *cos(w*t + sita);
            y_st(:,t+1) =  A0 *sin(w*t + sita);
        end
        %Ǳͧ��̽Ǳ����֮��ľ���
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
%���ָ���
A = find_cnt>0;
figure(1)
plot(x_t,y_t,x_st,y_st)
find_pro = 1*sum(A(:))/N
