%����delta-f��gk��fk�еĵ�Ч

clear;
close all;
%% ��������
%�������
m=1.67e-27; %kg, ��������
me=0.91e-30;
e=1.6e-19;  %C,���ӵ��
epsilon0=8.854e-12; %F, ��ռ�����
T0=150*e;%1500*e;   %J,�����¶ȳ��Բ�����������
v_th=sqrt(T0/m);    %�������ٶ�
n0=1e17;%1.0e19;  %1/m^3,�����ܶ�
%B0=1.0;     %T, �Ÿ�Ӧǿ��
%R0=1.0;     %m, ���ȵ�λ����Ϊ��GTC��λ�Ա�ʹ��
%phi0=1.0e-17*B0^2*R0^2*e/m; %V,����������ƵĴ�С
%omega=1.36911*e*B0/m;%1.54391*e*B0/m; %1/s����ӵĲ���Ƶ��
%x0=0.08660254;
%x1=0.1161895;

omegap=sqrt(n0*e^2/epsilon0/m);
omegape=sqrt(n0*e^2/epsilon0/me);
omegac=0.1*omegap;
B0=omegac*m/e;
omegace=e*B0/me;
krho=0.7;
k=krho/(m*v_th/e/B0);% �����Ĳ���
lambdaD=sqrt(epsilon0*T0/(n0)/e^2);
rho0=n0*e/epsilon0;
L=4*pi/k;
x0=0;
x1=L+x0;
%�������
xgrids=127;%һά�Ŀռ���
dx=L/xgrids;
xrange=x0:dx:x1; %�ռ䷶Χ
tstep=0.4/omegap;%0.5e-4*R0/v_th;%5.4143e-11;% %ʱ�䲽��
mstep=5*ceil(2*pi/omegac/tstep)
micell=200;
mi=xgrids*micell; % marker��Ŀ
deltan=zeros(1,xgrids); %�Ŷ��������ܶ�
deltaE=zeros(1,xgrids); %�Ŷ��糡
poissonE=zeros(xgrids,1); %��poisson���̽�����Ŷ��糡
nonlinear=1;
%�����
Etime=zeros(xgrids,mstep);
gyro_phi=zeros(xgrids,mstep);
ntime=zeros(xgrids,mstep);
parray=[1:mi 1:mi];
gk_phi_old = zeros(xgrids, mstep);
gk_phi_new = zeros(xgrids, mstep);
rhos = m*v_th/e/B0;
psir = @(r) (r.^2/2)/1000;
meshni=@(x) 0.7 - tanh(((psir(0.5*(x1+x0)-abs(x-0.5*(x1+x0))))/psiw- 0.15)*300)/4;

gx=rand(1,mi)'*(x1-x0)+x0;%linspace(x0,x1,mi)'; %���ǵ��ĵ�����
%gx=0.01*L*cos(k*gx)+gx;
vpx=v_th*randn(mi,1);
vpy=v_th*randn(mi,1);
pw=0.01*sin(k*gx);

rhox=m*vpy/e/B0;
xp=gx-rhox;

g0=floor(xp/dx-.5)+1;%floor((xp-x0)/dx-0.5)+1;
g=[g0;g0+1];
out=(g<1);g(out)=g(out)+xgrids;
out=(g>xgrids);g(out)=g(out)-xgrids;
h1=abs(xp/dx-g0+.5);%mod((xp-x0)/dx+0.5,1);
h=[1-h1;h1];
mat=sparse(parray,g,h,mi,xgrids);

dnf=(pw'*mat)*n0/micell;


rho_perp=m*sqrt(vpx.^2+vpy.^2)/e/B0;
%     gx=xp-rhox;
    gx_left=gx-rho_perp;
    gx_right=gx+rho_perp;
    
    g0=floor(gx/dx-.5)+1;%floor((xp-x0)/dx-0.5)+1;
%     g0_left=floor(gx_left/dx-.5)+1;
%     g0_right=floor(gx_right/dx-.5)+1;
    g=[g0;g0+1];
%     g_left=[g0_left;g0_left+1];
%     g_right=[g0_right;g0_right+1];
    out=(g<1);g(out)=g(out)+xgrids;
%     out_left=(g_left<1);g_left(out_left)=g_left(out_left)+xgrids;
%     out_right=(g_right<1);g_right(out_right)=g_right(out_right)+xgrids;
    out=(g>xgrids);g(out)=g(out)-xgrids;
%     out_left=(g_left>xgrids);g_left(out_left)=g_left(out_left)-xgrids;
%     out_right=(g_right>xgrids);g_right(out_right)=g_right(out_right)-xgrids;
    h1=abs(gx/dx-g0+.5);%mod((xp-x0)/dx+0.5,1);
%     h1_left=abs(gx_left/dx-g0_left+.5);
%     h1_right=abs(gx_right/dx-g0_right+.5);
    h=[1-h1;h1];
%     h_left=[1-h1_left;h1];
%     h_right=[1-h1_right;h1];
    mat=sparse(parray,g,h,mi,xgrids);
%     mat_left=sparse(parray,g_left,h_left,mi,xgrids);
%     mat_right=sparse(parray,g_right,h_right,mi,xgrids);
   
    
    g_rho= + ((pw.*besselj(0,k*rho_perp))'*mat)*n0/micell;%*0.5...
           %+0.25*((pw'*mat_left))*n0/micell...
           %+0.25*((pw'*mat_right))*n0/micell;%-rho0 * epsilon0;



figure;plot(dnf); hold on; plot(g_rho);%plot(0.01*n0*sin(k*xrange)/2.71828^2)