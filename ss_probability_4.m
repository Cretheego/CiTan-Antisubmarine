function  p=ss_probability_4(x,k,Vse,W,R0)
% t0 = s0/(Vse*sqrt(2/pi));
% t1 = D /(Vse + Vs);
% %t1 = D /(Vse + 300);
% R0 = Vse * (t1+t0); %ºŸ…Ëµƒ≥ı ºæ‡¿Î
u1 = Vse - W * Vse/(2*R0).*exp(-k*x);
if u1<0
   u1 = 0;
end
A = (u1>0);
u1 = u1.* A;
r1 = u1*
%u1 = u1.*exp(-k*x);
u2 = Vse + W * Vse/(2*R0)*exp(-k*x);
%p=1/(2*pi)*(pi*exp(-u1.^2)/(4*Vse^2)-(pi*exp(-u2.^2)/(4*Vse^2)));
p=1/(2*pi)*(pi*u2.*exp(-u2.^2)/(4*Vse^2)-(pi*u1.*exp(-u1.^2)/(4*Vse^2)));
end
