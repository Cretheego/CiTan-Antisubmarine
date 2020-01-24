clear all
%��������
%��ʼ��(m)
x0 = 0; y0 = 0;
%ɢ������
L = 4*10^3;
%���ִ���
%�������
h = 1;
N = 500/h;
find_cnt = zeros(1,N);
%��Ǳ�����ٶ�
V_s = 50*10^3;
%����Ƶ��
A0 = L;
w = 1;
%W_m�Ǵ�̽�ǵ��������
W_m = 600;
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
    x_st(:,1) = 2*L*cos(2*pi*(1/num:1/num:1))';
    y_st(:,1) = 2*L*sin(2*pi*(1/num:1/num:1))';   
    A0 = L;
    %�ٶ�(��) 5-13kn
    V = (5+(13-5)*rand(1))*2;
    %����[0,2*pi]���ȷֲ�
    H0 = 2*pi*rand(1);
    t = 0;
    sita = 0;
    %��С�ļ��
    A1 = (R_max - L) / T;
    while(1)
        x_t(t+1) = x_t(1) + V*cos(H0)*t;
        y_t(t+1) = y_t(1) + V*sin(H0)*t;
               
        %ʱ�̴�̽��t��λ��
        w = V_s / A0;

        x_st(:,t+1) = x_st(:,1) + A0 *cos(w*t + sita);
        y_st(:,t+1) = y_st(:,1) + A0 *sin(w*t + sita);

        A0 = A0 + 30;
        %Ǳͧ��̽Ǳ����֮��ľ���
        d_fq = sqrt((x_st(:,t+1) - repmat(x_t(t+1),num,1)).^2 + (y_st(:,t+1) - repmat(y_t(t+1),num,1)).^2);
        I = (d_fq <= W_m/2);
        if ~isempty(find(I(:)~=0,1))
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
plot(x_t,y_t,x_st(1,:),y_st(1,:))
find_pro = 1*sum(A(:))/N
