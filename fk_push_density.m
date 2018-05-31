%测试delta-f在gk和fk中的等效

clear;

%% 参数设置
%物理参数
m=1.67e-27; %kg, 粒子质量
me=0.91e-30;
e=1.6e-19;  %C,粒子电荷
epsilon0=8.854e-12; %F, 真空极化率
T0=150*e;%1500*e;   %J,粒子温度乘以玻尔兹曼常数
v_th=sqrt(T0/m);    %粒子热速度
n0=1e17;%1.0e19;  %1/m^3,粒子密度
%B0=1.0;     %T, 磁感应强度
%R0=1.0;     %m, 长度单位，仅为与GTC单位对比使用
%phi0=1.0e-17*B0^2*R0^2*e/m; %V,外加驱动电势的大小
%omega=1.36911*e*B0/m;%1.54391*e*B0/m; %1/s，外加的波的频率
%x0=0.08660254;
%x1=0.1161895;

omegap=sqrt(n0*e^2/epsilon0/m);
omegape=sqrt(n0*e^2/epsilon0/me);
omegac=0.1*omegap;
B0=omegac*m/e;
omegace=e*B0/me;
krho=0.43;
k=krho/(m*v_th/e/B0);% 波动的波数
lambdaD=sqrt(epsilon0*T0/(n0)/e^2);
rho0=n0*e/epsilon0;
L=4*pi/k;
x0=0;
x1=L+x0;
%程序参数
xgrids=127;%一维的空间格点
dx=L/xgrids;
xrange=x0:dx:x1; %空间范围
tstep=0.4/omegap;%0.5e-4*R0/v_th;%5.4143e-11;% %时间步长
mstep=5*ceil(2*pi/omegac/tstep);
micell=500;
mi=xgrids*micell; % marker数目
deltan=zeros(1,xgrids); %扰动的粒子密度
deltaE=zeros(1,xgrids); %扰动电场
poissonE=zeros(xgrids,1); %解poisson方程解出的扰动电场
nonlinear=1;
%诊断量
Etime=zeros(xgrids,mstep);
xtime=zeros(1,mstep);
vtime=zeros(1,mstep);
vytime=zeros(1,mstep);
gyro_phi=zeros(xgrids,mstep);
ntime=zeros(xgrids,mstep);
parray=[1:mi 1:mi];
gk_phi_old = zeros(xgrids, mstep);
gk_phi_new = zeros(xgrids, mstep);
rhos = m*v_th/e/B0;
psir = @(r) (r.^2/2)/1000;
psiw = 0.0375e-3;
meshn=@(x) ones(size(x));%0.7 - tanh(((psir(0.1183-abs(x/L*0.0634+0.0866-0.1183)))/psiw-0.15)*200)/10;
meshni=meshn(dx/2:dx:L-dx/2)';
meshti=ones(xgrids,1);
kapani =@(x) (meshn(x+0.001)-meshn(x-0.001))/0.001;
%矩阵
Mdr = zeros(xgrids);
for i=2:xgrids-1
    Mdr(i,i+1)=0.5/dx;
    Mdr(i,i-1)=-0.5/dx;
end

Mdr(1,end)=-0.5/dx;Mdr(1,2)=0.5/dx;
Mdr(end,1)=0.5/dx; Mdr(end, end-1) = -0.5/dx;

Mlap=zeros(xgrids);
for i=2:xgrids-1
    Mlap(i,i)=-2/dx^2;
    Mlap(i,i+1)=1/dx^2;
    Mlap(i,i-1)=1/dx^2;
end
Mlap(1,1)=-2/dx^2;
Mlap(1,2)=1/dx^2;
Mlap(1,xgrids)=1/dx^2;
Mlap(end,end)=-2/dx^2;
Mlap(end,end-1) = 1/dx^2;
Mlap(end,1) = 1/dx^2;

Mlap_new = Mlap;
%Mlap(end,:) = 1;

for i=1:xgrids
    Mlap_new(i,:) = Mlap_new(i,:)+kapani(dx/2+i*dx)*Mdr(i,:);
end

Mold = Mlap;
for i=1:xgrids-1
    Mold(i,:) = -Mold(i,:)*e*(n0*meshni(i))/(T0*meshti(i))*rhos^2;
end
Mold = (eye(xgrids)-rhos^2*Mlap)\Mold;%-Mlap;
Mold(end,:)=1;
% Mold(1,:)=0;
% Mold(1,1)=1;
% Mold(end,:)=0;
% Mold(end,end)=1;

%Mlap_new=Mlap;
Mnew = Mlap_new;
for i=1:xgrids-1
    Mnew(i,:) = -Mnew(i,:)*e*(n0*meshni(i))/(T0*meshti(i))*rhos^2;
end
Mnew = (eye(xgrids)- rhos^2*Mlap_new)\Mnew;%-Mlap;
Mnew(end,:)=1;

%开始
gx=linspace(x0,x1,mi)';%rand(1,mi)'*(x1-x0)+x0;% %这是导心的坐标
%gx=0.01*L*cos(k*gx)+gx;
vpx=v_th*randn(mi,1);
vpy=v_th*randn(mi,1);


rhox=m*vpy/e/B0;
xp=gx-rhox;
%gx=xp+rhox;
gw=0.0*sin(k*gx);

g0=floor(xp/dx-.5)+1;%floor((xp-x0)/dx-0.5)+1;
g=[g0;g0+1];
out=(g<1);g(out)=g(out)+xgrids;
out=(g>xgrids);g(out)=g(out)-xgrids;
h1=abs(xp/dx-g0+.5);%mod((xp-x0)/dx+0.5,1);
h=[1-h1;h1];
mat=sparse(parray,g,h,mi,xgrids);
pw=gw;
dnf=(pw'*mat)*n0/micell;


rho_perp=m*sqrt(vpx.^2+vpy.^2)/e/B0;

g_g0=floor(gx/dx-.5)+1;
g_g=[g_g0;g_g0+1];
out=(g_g<1);g_g(out)=g_g(out)+xgrids;
out=(g_g>xgrids);g_g(out)=g_g(out)-xgrids;
g_h1=abs(gx/dx-g_g0+.5);
g_h=[1-g_h1;g_h1];
g_mat=sparse(parray,g_g,g_h,mi,xgrids);
g_rho= + ((gw.*besselj(0,k*rho_perp))'*g_mat)*n0/micell;%*0.5...
%+0.25*((pw'*mat_left))*n0/micell...
%+0.25*((pw'*mat_right))*n0/micell;%-rho0 * epsilon0;

dtime=tstep;
figure;
extphi=1.0e4*sin(k*[dx/2:dx:L-dx/2]');
for istep=1:mstep
    t=dtime*istep;
    extphi=1.0e5*sin(k*[dx/2:dx:L-dx/2]')*sin(0.06*istep*omegac*tstep);
    extE=-Mdr*extphi;
    xp=xp./L+10.0; xp=L.*(xp-floor(xp));
    %push x
    xp=xp+vpx*dtime;
%     xp=gx-rhox;
    xtime(istep)=xp(15000);
    vtime(istep)=vpx(15000);
    vytime(istep)=vpy(15000);
    xp=xp./L+10.0; xp=L.*(xp-floor(xp));
    
    g0=floor(xp/dx-.5)+1;%floor((xp-x0)/dx-0.5)+1;
    g=[g0;g0+1];
    out=(g<1);g(out)=g(out)+xgrids;
    out=(g>xgrids);g(out)=g(out)-xgrids;
    h1=abs(xp/dx-g0+.5);%mod((xp-x0)/dx+0.5,1); 
    h=[1-h1;h1];
    mat=sparse(parray,g,h,mi,xgrids);

    dnf=((meshn(gx).*pw)'*mat)*n0/micell;
    pw=pw+vpx.*(mat*extE)*e/T0*dtime;
    gammab=e*B0/m*dtime/(1+e^2*B0^2/m^2/4*dtime^2);
    vx1=vpx+(vpy-vpx*e*B0/m*dtime/2)*gammab;
    vpy=vpy-(vpx+vpy*e*B0/m*dtime/2)*gammab;
    vpx=vx1;
    
    rhox=m*vpy/e/B0;
    gx=xp+rhox;
    g_g0=floor(gx/dx-.5)+1;
    g_g=[g_g0;g_g0+1];
    out=(g_g<1);g_g(out)=g_g(out)+xgrids;
    out=(g_g>xgrids);g_g(out)=g_g(out)-xgrids;
    g_h1=abs(gx/dx-g_g0+.5);
    g_h=[1-g_h1;g_h1];
    g_mat=sparse(parray,g_g,g_h,mi,xgrids);
%     g_rho= + ((gw.*besselj(0,k*rho_perp).*meshni(gx))'*g_mat)*n0/micell;
    g_rho= + ((pw.*meshn(gx))'*g_mat)*n0/micell;
    
    
    %guiding center dwgdt=0; dxdt=0;
    subplot(221);plot(extphi);
    subplot(222);plot(dx/2:dx:L-dx/2,dnf,dx/2:dx:L-dx/2,g_rho);
    drawnow;
end
figure;plot(g_rho);hold on; plot(dnf);
figure;plot(g_rho-dnf); hold on;
plot(Mold*extphi);
plot(Mnew*extphi);
% fix=e/T0*n0.*meshni*krho^2/k.*kapani(dx/2:dx:L-dx/2)'*(besselj(0,krho^2)-besselj(1,krho^2))*exp(-krho^2)*1.0e5;
% plot(e*extphi/T0*n0.*meshni*(1-besselj(0,krho^2)*exp(-krho^2)));
% plot(e*extphi/T0*n0.*meshni*(1-besselj(0,krho^2)*exp(-krho^2))-(fix.*cos(k*[dx/2:dx:L-dx/2]')));
%plot(0.01*n0*sin(k*xrange)/2.71828^2)