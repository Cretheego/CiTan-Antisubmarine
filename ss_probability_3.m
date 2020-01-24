function  p=ss_probability_3(r,x,ss,k,Vse,Vs,W,D,s0)
D1 = sqrt(D^2+r.^2-2*r.*D.*cos(ss));
t1 = D1 ./(Vse + Vs);
R0 = Vse * t1; %ºŸ…Ëµƒ≥ı ºæ‡¿Î
u1 = Vse - W * Vse./(2*R0).*exp(-k*x);
if u1<0
   u1 = 0;
end
A = (u1>0);
u1 = u1.* A;
%u1 = u1.*exp(-k*x);
u2 = Vse + W * Vse./(2*R0).*exp(-k*x);
p=1/(2*pi)^2*(pi*exp(-u1.^2)/(4*Vse^2)-(pi*exp(-u2.^2)/(4*Vse^2))).*(r/s0^2).*exp(-r.^2/(2*s0^2));
end
