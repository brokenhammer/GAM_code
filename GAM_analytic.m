%使用Gao 2006 文章中的解析表达式求解色散关系
%% physical parameters
mi=1;
me=1/1836;
vti=1;
tau=1;
q=1.4;
ep=0.2;
B0=1;
e=1;
rhoti=mi*vti/e/B0;
omega0=e*B0/mi;
b=0.2^2;
kr=sqrt(b)/rhoti;
Ti=mi*vti^2;
Te=Ti*tau;
R0=1;
n0=1;
syms vpa vpe theta omega;
j0=besselj(0, kr*vpe/omega0);
omegat = vpa/q/R0;
Fm = n0/sqrt(2*pi)^3/vti^3*exp(-vpa^2/2/vti^2)*exp(-vpe^2/2/vti^2);
omegad = kr/R0/omega0*(1/2*vpe^2+vpa^2);
g = 0;
for n=-6:6
    g = g+e*Fm/Ti*j0*(omega)*besselj(n,omegad/omegat)^2/(omega+n*omegat);
end
g = g - e*Fm/Ti*j0*(omega);

solo = solve(e*n0/Ti*(1-besselj(0,b)*exp(-b))==2*pi*int(vpe*int(g,vpa,-inf,+inf),vpe,0,inf),omega);
