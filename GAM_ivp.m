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
vpara_grid = 99;
vperp_grid = 100;
theta_grid = 101;
vpara_array = vti*linspace(-5,5,vpara_grid);
dvpa = vpara_array(2) - vpara_array(1);
vperp_array = vti*linspace(0,5,vperp_grid);
dvpe = vperp_array(2) - vperp_array(1);
theta_array = linspace(-10,10,theta_grid);
dtheta = theta_array(2) - theta_array(1);
niter = 2;
[vpara,tspace,vperp] = meshgrid(vpara_array, theta_array, vperp_array);

dt = 0.1;
nt=20;
%% initialize
g0 = 0.01 * exp(-vpara.^2/2/vti^2) .* exp(-vperp.^2/2/vti^2) .* exp(-(tspace/2).^2) + 0 * 1i;
omega_di = -kr*rhoti/vti/R0*(0.5*vperp.^2+vpara.^2).*sin(tspace);
Fm = n0/(sqrt(2*pi))^3/vti^3*exp(-vpara.^2/2/vti^2) .* exp(-vperp.^2/2/vti^2);
gamma0 = besseli(0,kr^2*rhoti^2)*exp(-kr^2*rhoti^2);
Mdt = zeros(theta_grid);
for k=2:theta_grid-1
    Mdt(k,k-1) = -1/2/dtheta;
    Mdt(k,k+1) = 1/2/dtheta; 
end
Mdt(1) = 1;
Mdt(end) = 1;

% Mpos = zeros(vpara_grid, vperp_grid, theta_grid, theta_grid);
% Mneg = zeros(vpara_grid, vperp_grid, theta_grid, theta_grid);
% Mpush = zeros(vpara_grid, vperp_grid, theta_grid, theta_grid);
% % iterative
% % Mpos * g+ = Mneg * g- + (-Mphi-Mgphi)*phi(1.5)
% % Mpoisson * phi(1.5) = int{J0*(g+ + g-)/2}
% tic;
% for ivp = 1:vpara_grid
%     for ivperp = 1:vperp_grid
%         Mpartg = Mdt*vpara_array(ivp)/q./(1.0+ep*cos(theta_array))';
%         Mg = diag(1/dt + 1i/2*omega_di(:,ivp,ivperp));
%         Mpos(ivp, ivperp,:,:) = Mpartg + Mg;
%         Mg = diag(1/dt - 1i/2*omega_di(:,ivp,ivperp));
%         Mneg(ivp, ivperp,:,:) = -Mpartg + Mg;
%         Mpartphi = Mpartg*e/Ti*besselj(0,kr*vperp_array(ivperp)/omega0).*Fm(:,ivp,ivperp);
%         Mphi = diag(1i*omega_di(:,ivp,ivperp)*e/Ti*besselj(0,kr*vperp_array(ivperp)/omega0).*Fm(:,ivp,ivperp));
%         Mpush(ivp,ivperp,:,:) = -Mpartphi-Mphi;
%     end
%     if(mod(ivp,10)==0)
%         ivp
%     end
% end
% Mpos(:,:,1,:) = 0;
% Mpos(:,:,1,1) = 1;
% Mpos(:,:,end,:) = 0;
% Mpos(:,:,end,end) = 1;
% toc;
Mpoisson = e*n0/Te*(eye(theta_grid))+e*n0/Ti*(eye(theta_grid)*(1-gamma0));

% onestep
% tic;

g_intvpa = squeeze(sum(g0, 2))*dvpa;
g_intvpe = sum(g_intvpa.*vperp_array.*besselj(0,kr*vperp_array/omega0), 2)*dvpe;
deltan0 = 2*pi*g_intvpe;
phi0 = Mpoisson\deltan0;

phi_track = zeros(theta_grid,nt);
g_neg = g0;
g_pos = g0;
it=1;
for t = 0:dt:dt*(nt-1)
for iter = 1:niter
    g_mid = (g_neg+g_pos)/2;
    g_intvpa = squeeze(sum(g_mid, 2))*dvpa;
    g_intvpe = sum(g_intvpa.*vperp_array.*besselj(0,kr*vperp_array/omega0), 2)*dvpe;
    deltan = 2*pi*g_intvpe;
    phi_mid = Mpoisson\deltan;
    for ivp = 1:vpara_grid
        for ivperp = 1:vperp_grid
            g_pushed = (Mneg(ivp,ivperp)*g_neg(:,ivp,ivperp)+Mpush(ivp,ivperp)*phi_mid);
            g_pushed(1) = 0;
            g_pushed(end) = 0;
            g_pos(:,ivp,ivperp) = Mpos(ivp,ivperp)\g_pushed;
        end
    end
end
phi_track(:,it)=phi_mid;
it=it+1;
g_neg = g_pos;
end
% toc;

figure;plot(theta_array,phi0);
figure;subplot(211);plot(theta_array,real(phi_mid));subplot(212);plot(theta_array,imag(phi_mid));
figure;plot(abs(phi_track(ceil(end/2),:)));