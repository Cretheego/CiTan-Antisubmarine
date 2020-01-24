function  p=ss_probability(x,k,Vse,W,R0)
% t1 = D /(Vse + Vs);
% R0 = Vse * t1;
% k = tan(asin(Vse/Vs));
u1 = Vse - W * Vse/(2*R0).*exp(-k*x);
if u1<0
   u1 = 0;
end
A = (u1>0);
u1 = u1.* A;
%u1 = u1.*exp(-k*x);
u2 = Vse + W * Vse/(2*R0)*exp(-k*x);
p=1/(2*pi)*(pi*exp(-u1.^2)/(4*Vse^2)-(pi*exp(-u2.^2)/(4*Vse^2)));
end
