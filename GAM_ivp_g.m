% 1��ʹ��initial value problem�ķ������GAMƵ�ʺ������ʣ����ظ������������ĵ�����
% 2��˳����֤Gao 2006�еĽ�����
% 3����֤Te->0 limit,������ƽ������
% 4��֮������omega(*T),ʹ����ֵ���ͽ������������������ITG���
% 5����Poisson��������������
% 6��ʹ������������ģ�任���̣�����ƽ���任����ʲô��

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
g0 = 0.01 * exp(-vpara.^2/2/vti^2) .* exp(-vperp.^2/2/vti^2) .* exp(-(tspace/2).^2) + 0 * 1i;
% h0 = 0.01 * exp(-vpara.^2/2/vti^2) .* exp(-vperp.^2/2/vti^2) .* exp(-(tspace/2).^2) + 0 * 1i;
omega_di = -kr*rhoti/vti/R0*(0.5*vperp.^2+vpara.^2).*sin(tspace);
Fm = n0/(sqrt(2*pi))^3/vti^3*exp(-vpara.^2/2/vti^2) .* exp(-vperp.^2/2/vti^2);
gamma0 = besseli(0,kr^2*rhoti^2)*exp(-kr^2*rhoti^2);
deltan = zeros(theta_grid);
pgpth = zeros(theta_grid, vpara_grid, vperp_grid);
dgdt1 = pgpth;
dgdt2 = pgpth;
%RK 4 for h=h+dhdt*dt

for it = 1:nt
    g = g0;
    %calc dgdt1
    % integrate for phi
    phi = calc_phi_g(g,vperp,omega0,dvpe,dvpa,kr);
    dgdt1 = gtime(g,phi,omega_di,q,R0,pgpth,kr,vperp,omega0,vpara);
    
    %calc dgdt2
    g = g0 + dgdt1*dt/2;
    