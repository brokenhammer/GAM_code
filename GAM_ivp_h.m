% 1、使用initial value problem的方法求解GAM频率和增长率，先重复出仇智勇论文的内容
% 2、顺便验证Gao 2006中的解析解
% 3、验证Te->0 limit,即磁面平均方程
% 4、之后增加omega(*T),使用数值（和解析方法？）求解增加ITG后的
% 5、在Poisson方程增加修正项
% 6、使用完整的气球模变换方程？磁面平均变换后是什么？

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
kr=0.2/rhoti;
Ti=mi*vti^2;
Te=Ti*tau;
R0=1;
n0=1;

%% control parameters
vpara_grid = 21;
vperp_grid = 21;
theta_grid = 121;
vpara_array = vti*linspace(-5,5,vpara_grid);
dvpa = vpara_array(2) - vpara_array(1);
vperp_array = vti*linspace(0,5,vperp_grid);
dvpe = vperp_array(2) - vperp_array(1);
theta_array = linspace(-2*pi,2*pi,theta_grid);
dtheta = theta_array(2) - theta_array(1);
niter = 2;
[vpara,tspace,vperp] = meshgrid(vpara_array, theta_array, vperp_array);

dt = 0.1;
nt=20;
%% initialize
%g0 = 0.01 * exp(-vpara.^2/2/vti^2) .* exp(-vperp.^2/2/vti^2) .* exp(-(tspace/2).^2) + 0 * 1i;
h0 = 0.01 * exp(-vpara.^2/2/vti^2) .* exp(-vperp.^2/2/vti^2) .* exp(-(tspace/2).^2) + 0 * 1i;
omega_di = -kr*rhoti/vti/R0*(0.5*vperp.^2+vpara.^2).*sin(tspace);
Fm = n0/(sqrt(2*pi))^3/vti^3*exp(-vpara.^2/2/vti^2) .* exp(-vperp.^2/2/vti^2);
gamma0 = besseli(0,kr^2*rhoti^2)*exp(-kr^2*rhoti^2);
deltan = zeros(theta_grid);
%RK 4 for h=h+dhdt*dt
for it = 1:nt
    % integrate for phi
    for ith = theta_grid
    for ivp = 1:vpara_grid
        for ivc = 1:vperp_grid
            vpe = vperp_array(ivc);
            krrho = kr * vpe / omega0;
            deltan(ith) = dletan(ith) + vpe * h0(ith,icp,ivc) * dvpe * dvpa * besselj(0,krrho);
        end
    end
    end
    
    phi = deltan ./ (e * n0 / Ti);
    
    