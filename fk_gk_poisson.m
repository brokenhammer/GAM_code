%1d-3v的delta-f，静电code，试图解出带回旋运动的色散关系（国际单位制）
clear;
close all;
%% 参数设置
%物理参数
m=1.67e-27; %kg, 粒子质量
me=0.91e-30;
e=1.6e-19;  %C,粒子电荷
epsilon0=8.854e-12; %F, 真空极化率
T0=150*e;%1500*e;   %J,粒子温度乘以玻尔兹曼常数
v_th=sqrt(T0/m);    %粒子热速度
n0=1e18;%1.0e19;  %1/m^3,粒子密度
%B0=1.0;     %T, 磁感应强度
%R0=1.0;     %m, 长度单位，仅为与GTC单位对比使用
%phi0=1.0e-17*B0^2*R0^2*e/m; %V,外加驱动电势的大小
%omega=1.36911*e*B0/m;%1.54391*e*B0/m; %1/s，外加的波的频率
%x0=0.08660254;
%x1=0.1161895;

omegap=sqrt(n0*e^2/epsilon0/m);
omegape=sqrt(n0*e^2/epsilon0/me);
omegac=0.6*omegap;
B0=omegac*m/e;
omegace=e*B0/me;
k=0.2/(m*v_th/e/B0);% 波动的波数
lambdaD=sqrt(epsilon0*T0/(n0)/e^2);
rho0=n0*e/epsilon0;
L=2*pi*2/k;
x0=0;
x1=L+x0;
%程序参数
xgrids=101;%一维的空间格点
dx=L/xgrids;
xrange=x0:dx:x1; %空间范围
tstep=0.1/omegap;%0.5e-4*R0/v_th;%5.4143e-11;% %时间步长
mstep=5*ceil(2*pi/omegac/tstep)
micell=200;
% zion=zeros(5,xgrids*micell); %粒子信息 1维位置，3个速度分量，1个权重
% vgrids=201;
% vrange=[-7*v_th,7*v_th,vgrids];
% dv=vrange(2)-vrange(1);
mi=xgrids*micell; % marker数目
deltan=zeros(1,xgrids); %扰动的粒子密度
deltaE=zeros(1,xgrids); %扰动电场
poissonE=zeros(xgrids,1); %解poisson方程解出的扰动电场
nonlinear=1;
%诊断量
Etime=zeros(xgrids,mstep);
gyro_phi=zeros(xgrids,mstep);
ntime=zeros(xgrids,mstep);
parray=[1:mi 1:mi];
gk_phi_old = zeros(xgrids, mstep);
gk_phi_new = zeros(xgrids, mstep);
rhos = m*v_th/e/B0;

%% 粒子初始化
%rng('default'); %随机数种子设为default，每次运行结果相同，注释掉则以时间为种子
xp=linspace(x0,x1,mi)';%rand(1,mi)'*(x1-x0)+x0; %粒子位置，均匀随机分布%
% xp=0.01*L*cos(k*xp)+xp; %初始扰动
vpx=v_th*randn(mi,1);%normrnd(0.0,v_th,[mi,1]);%粒子速度，玻尔兹曼分布,注意与GTC的差别？
vpy=v_th*randn(mi,1);%normrnd(0.0,v_th,[mi,1]);
vpz=v_th*randn(mi,1);%normrnd(0.0,v_th,[mi,1]);
pw=zeros(mi,1)+1;

A=zeros(xgrids,xgrids);
for i=2:xgrids-1
    A(i,i)=2/dx^2;
    A(i,i-1)=-1/dx^2;
    A(i,i+1)=-1/dx^2;
end
A(1,1)=2/dx^2;
A(1,2)=-1/dx^2;
A(1,xgrids)=-1/dx^2;
B=A;
A(end,:)=1;
A=A*n0*(m/B0^2);
B(end,1)=-1/dx^2;
B(end,end)=2/dx^2;
B(end,end-1)=-1/dx^2;
C=B;
B=B*lambdaD^2;
B(end,:)=1;
Diff=zeros(xgrids);
for i=2:xgrids-1
    Diff(i,i+1)=1/2/dx;
    Diff(i,i-1)=-1/2/dx;
end
Diff(end,:) = 1;
Diff(1,1) = 1;
% Diff(1,2)=1/2/dx;
% Diff(1,end)=-1/2/dx;
% Diff(end,1)=1/2/dx;
% Diff(end,end-1)=-1/2/dx;


a=1;
psir = @(r) (r.^2/2)/1000;
psiw = 0.0375e-3;
psi0 = 0.1*psiw;
psi1 = 0.3*psiw;
rpsi = @(psi) sqrt(2*psi*1000);
rmesh = linspace(rpsi(psi0), rpsi(psi1), xgrids)';
ppsi = psir(xp / L * (rpsi(psi1) - rpsi(psi0)) + rpsi(psi0));
for i = 1:mi
    if xp(i) > 0.5 * L
        ppsi(i) = psir((1-xp(i) / L) * (rpsi(psi1) - rpsi(psi0)) + rpsi(psi0));
    end
end
psimesh=psir(rmesh);
meshne = zeros(xgrids,1);
meshne(1:ceil(xgrids/2)) = 1.0;%0.7 - tanh(((psimesh(1:ceil(xgrids/2))/psiw) - 0.15)*200)/10;
meshne(ceil(xgrids/2):xgrids) = meshne(ceil(xgrids/2):-1:1);
figure;plot(meshne);

% for i=1:mi
%     pw(i) = pw(i) * (0.7 - tanh((ppsi(i)/psiw - 0.12)*45)/10)/meshne(1);
% end
pw = pw+0.1*cos(k*xp);
meshne = meshne / meshne(1);
meshni=meshne;
meshte=ones(xgrids,1);
meshti=meshte;
% rho0= rho0 * meshni';
%dr = rmesh(1) - rmesh(0);
Mdr = zeros(xgrids);
for i=2:xgrids-1
    Mdr(i,i+1)=0.5/dx;
    Mdr(i,i-1)=-0.5/dx;
end

Mdr(1,end)=-0.5/dx;Mdr(1,2)=0.5/dx;
Mdr(end,1)=0.5/dx; Mdr(end, end-1) = -0.5/dx;

kapani=Mdr*meshni./meshni;


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
    Mlap_new(i,:) = Mlap_new(i,:)+1*kapani(i)*Mdr(i,:);
end
% Mlap_new(1,1)=-2/dx^2;
% Mlap_new(1,2)=1/dx^2;
% Mlap_new(1,xgrids)=1/dx^2;
% Mlap_new(end,:)=1;

% Mold = Mlap;
% for i=1:xgrids -1
%     Mold(i,:) = -Mold(i,:)*(e^2*(n0*meshni(i))/(meshti(i)*T0/rhos^2));
% end

Mold = Mlap;
for i=1:xgrids-1
    Mold(i,:) = -Mold(i,:)*epsilon0*meshti(i)*T0/e^2/(n0*meshni(i));
end
Mold = -(rhos^2*Mlap)+(eye(xgrids)-rhos^2*Mlap)*Mold;
% Mold(1,1)=-2/dx^2;
% Mold(1,2)=1/dx^2;
% Mold(1,xgrids)=1/dx^2;
Mold(end,:)=1;


% Mnew = Mlap_new;
% for i=1:xgrids-1
%     Mnew(i,:) = -Mnew(i,:)*(e^2*(n0*meshni(i))/(meshti(i)*T0/rhos^2));
% end
Mnew = Mlap;
for i=1:xgrids-1
    Mnew(i,:) = -Mnew(i,:)*epsilon0*meshti(i)*T0/e^2/(n0*meshni(i));
end

Mnew = -(rhos^2*Mlap_new)+(eye(xgrids)-rhos^2*Mlap_new)*Mnew;

% Mnew(1,1)=2/dx^2;
% Mnew(1,2)=-1/dx^2;
% Mnew(1,xgrids)=-1/dx^2;
Mnew(end,:)=1;

B=Mlap;
B(end,:)=1;


lambdaDe=lambdaD;
C=C/(lambdaDe^2+lambdaD^2);
for i=1:xgrids
    C(i,i)=C(i,i)+1/lambdaD^2;
end
C(end,:)=1;

%%
figure;
dtime=tstep;
for istep=1:mstep;
%     if(istep<=mstep)
%         extE=5.0e4*sin(k*[dx/2:dx:L-dx/2]-1.92*omegac*tstep*istep);
%     else
        extE=0.0;
%     end
%     extE=0;%phi0*k*sin(k*xrange-omega*istep*tstep);
%     phileft=0;
%     phiright=0;

%     deltaE=extE+poissonE;

    xp=xp./L+10.0; xp=L.*(xp-floor(xp));
    %push x
    xp=xp+vpx*dtime;
%     xp=mod(xp-x0,(x1-x0))+x0;%更新坐标，使用了周期性边界条件
    
    g0=floor(xp/dx-.5)+1;%floor((xp-x0)/dx-0.5)+1;
    g=[g0;g0+1];
    out=(g<1);g(out)=g(out)+xgrids;
    out=(g>xgrids);g(out)=g(out)-xgrids;
    h1=abs(xp/dx-g0+.5);%mod((xp-x0)/dx+0.5,1); 
    h=[1-h1;h1];
    mat=sparse(parray,g,h,mi,xgrids);

    rho=full(e/epsilon0*(pw'*mat))*n0/micell-rho0;
%     rho=full(e/epsilon0*(pw'*mat))*n0/micell;
%     pw = pw-mean(pw);
%     rho=0.05*rho0*cos(k*[dx/2:dx:L-dx/2]);

%     phi = rho'/k^2;
%     poissonE = -Mdr *phi;
%     poissonE=poissonE+extE';
    rho(end) = 0;
    %rho(1) = 0;
    %phi=B\(lambdaD^2*(rho'.*meshni));
    phi = -B\(rho' .* meshni);
    poissonE=-Mdr*phi;


    vpx=vpx+mat*poissonE*e/m*dtime/2;
    gammab=e*B0/m*dtime/(1+e^2*B0^2/m^2/4*dtime^2);
    vx1=vpx+(vpy-vpx*e*B0/m*dtime/2)*gammab;
    vpy=vpy-(vpx+vpy*e*B0/m*dtime/2)*gammab;
    vpx=vx1;
    vpx=vx1+mat*poissonE*e/m*dtime/2;
    

%     vpx=vpx+e/m*(mat*poissonE-vpz*B0)*dtime;
%     vpy=vpy-e/m*(vpx*B0)*dtime;
%     vpz=vpz+e/m*(vpx*B0)*dtime;
%     pw=pw+e/T0*mat*poissonE.*vpx*dtime;
    
    
    %rho_y=-m*vpx/e/B0;
    rho_x=m*vpy/e/B0;
    rho_perp=m*sqrt(vpx.^2+vpy.^2)/e/B0;
    gx=xp+rho_x;
    gx_left=gx-rho_perp;
    gx_right=gx+rho_perp;
    
    g0=floor(gx/dx-.5)+1;%floor((xp-x0)/dx-0.5)+1;
    g0_left=floor(gx_left/dx-.5)+1;
    g0_right=floor(gx_right/dx-.5)+1;
    g=[g0;g0+1];
    g_left=[g0_left;g0_left+1];
    g_right=[g0_right;g0_right+1];
    out=(g<1);g(out)=g(out)+xgrids;
    out_left=(g_left<1);g_left(out_left)=g_left(out_left)+xgrids;
    out_right=(g_right<1);g_right(out_right)=g_right(out_right)+xgrids;
    out=(g>xgrids);g(out)=g(out)-xgrids;
    out_left=(g_left>xgrids);g_left(out_left)=g_left(out_left)-xgrids;
    out_right=(g_right>xgrids);g_right(out_right)=g_right(out_right)-xgrids;
    h1=abs(gx/dx-g0+.5);%mod((xp-x0)/dx+0.5,1);
    h1_left=abs(gx_left/dx-g0_left+.5);
    h1_right=abs(gx_right/dx-g0_right+.5);
    h=[1-h1;h1];
    h_left=[1-h1_left;h1];
    h_right=[1-h1_right;h1];
    mat=sparse(parray,g,h,mi,xgrids);
    mat_left=sparse(parray,g_left,h_left,mi,xgrids);
    mat_right=sparse(parray,g_right,h_right,mi,xgrids);
    
%     g_rho = n0/micell*(full(e/epsilon0*(pw'*mat))*epsilon0 *0.5...
%         +0.25*full(e/epsilon0*(pw'*mat_left))*epsilon0...
%         +0.25*full(e/epsilon0*(pw'*mat_right))*epsilon0);
    
    g_rho=-rho0 * epsilon0 + full(e/epsilon0*(pw'*mat))*n0/micell*epsilon0;%*0.5...
          %+0.25*full(e/epsilon0*(pw'*mat_left))*n0/micell*epsilon0...
          %+0.25*full(e/epsilon0*(pw'*mat_right))*n0/micell*epsilon0;
    %g_rho=g_rho-mean(g_rho);

    %g_rho=0.05*rho0*epsilon0*cos(k*[dx/2:dx:L-dx/2]);
    g_rho1=g_rho;
    %g_rho1(1) = 0;
    g_rho1(end)=0;
    %g_phi=A\g_rho1';
    %g_phi=g_phi;%+g_rho'/(n0*e^2/T0);
    %g_phi=g_phi/(1+omegape^2/omegace^2);
    
%     rho_old = (eye(xgrids) + Mlap*rhos^2) * (g_rho'.*meshni);
    rho_old = (eye(xgrids) - Mlap*rhos^2)*meshti*T0/e^2/(n0).*(g_rho');
    rho_old(end) = 0;
    phi_old = Mold\rho_old;
    %phi_old = Mold\(g_rho1' .*meshni)+(g_rho1' .*meshni) /(n0*e^2/T0);%- 5.33* (g_rho1' .*meshni) /(n0*e^2/T0);
    
    %rho_new = (eye(xgrids) + Mlap_new*rhos^2) * (g_rho'.*meshni);
    rho_new = (eye(xgrids)-Mlap_new*rhos^2)*meshti*T0/e^2/(n0).*(g_rho');
    %rho_new = (rho_new + rho_new(end:-1:1))/2;

    rho_new(end) = 0;
    phi_new = Mnew\(rho_new);
    %phi_new = Mnew\(g_rho1' .*meshni)+(g_rho1' .*meshni) /(n0*e^2/T0);%- 5.33* (g_rho1' .*meshni) /(n0*e^2/T0);
    
    %diagnosis
%     if istep>mstep
    Etime(:,istep)=phi;
    %gyro_phi(:,istep)=g_phi;
%     Etime(:,istep)=poissonE;
    ntime(:,istep)=rho;
    gk_phi_old(:,istep)=phi_old;
    gk_phi_new(:,istep)=phi_new;
%     end

    subplot(3,2,1);
    plot(dx/2:dx:L-dx/2,poissonE);
    title(['t=',num2str(istep)]);
    %ylim([-1e6,1e6]);
    subplot(3,2,3);plot(dx/2:dx:L-dx/2,rho*epsilon0/e/n0);
    subplot(3,2,2);plot(dx/2:dx:L-dx/2,g_rho/e/n0);
    subplot(3,2,5);plot(dx/2:dx:L-dx/2,phi);
    subplot(3,2,4);plot(dx/2:dx:L-dx/2,phi_old');
    subplot(3,2,6);plot(dx/2:dx:L-dx/2,phi_new');
    drawnow;
end


% figure;
% fE=fft(Etime(xgrids/2,:));
% ft=[1:mstep]/mstep*(2*pi/(omegac*dtime));
% plot(ft,abs(fE));
% xlim([0,10])
figure;
mean_E=mean(Etime,2);
mean_g=mean(gyro_phi,2);
mean_gk_old = mean(gk_phi_old,2);
mean_gk_new = mean(gk_phi_new,2);
mean_n=mean(ntime,2);
subplot(211);
plot(mean_E,'k');
hold on;
plot(phi_old,'r--'); plot(phi_new,'b');hold off;
subplot(212);
plot(meshni);
%plot(Etime(xgrids/2,:));

