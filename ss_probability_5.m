function  P=ss_probability_5(k,Vse,W,R0,t0,f,N)

% A = (u1>0);
% u1 = u1.* A;
% %u1 = u1.*exp(-k*x);
% 
% %p=1/(2*pi)*(pi*exp(-u1.^2)/(4*Vse^2)-(pi*exp(-u2.^2)/(4*Vse^2)));
% %p=1/(4*Vse^2)*(pi*u2.*exp(-u2.^2/(4*Vse^2/pi))-(pi*u1.*exp(-u1.^2/(4*Vse^2/pi))));
% p=integral(@(y) 1/(4*Vse^2)*(pi*y.*exp(-y.^2/(4*Vse^2/pi))),u1,u2);
% end
syms x y u1 u2
u1 = Vse - W * Vse/(2*R0).*exp(-k*x);
u1 = Vse - W * Vse/(2*R0).*exp(-k*x);
u1 = (u1+abs(u1))/2;
u2 = Vse + W * Vse/(2*R0)*exp(-k*x);
sigma0 = 1000^2;
sigma2 = sigma0/t0^2+ Vse^2*2/pi;
%fun = 1/(4*Vse^2)*(y*exp(-y^2/(4*Vse^2/pi)));
fun = 1/(2*pi*sigma2)*(y*exp(-y^2/(2*sigma2)));
if f==0
    P = int(int(fun,y,u1,u2),x,0,2*pi);
else
    P = int(int(fun,y,u1,u2),x,0,2*N*pi);
end

