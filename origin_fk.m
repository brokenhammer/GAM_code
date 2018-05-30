%1d-3v��delta-f������code����ͼ����������˶���ɫɢ��ϵ�����ʵ�λ�ƣ�
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
n0=1e18;%1.0e19;  %1/m^3,�����ܶ�
%B0=1.0;     %T, �Ÿ�Ӧǿ��
%R0=1.0;     %m, ���ȵ�λ����Ϊ��GTC��λ�Ա�ʹ��
%phi0=1.0e-17*B0^2*R0^2*e/m; %V,����������ƵĴ�С
%omega=1.36911*e*B0/m;%1.54391*e*B0/m; %1/s����ӵĲ���Ƶ��
%x0=0.08660254;
%x1=0.1161895;
rho0=n0*e/epsilon0;
omegap=sqrt(n0*e^2/epsilon0/m);
omegape=sqrt(n0*e^2/epsilon0/me);
omegac=0.4*omegap;
B0=omegac*m/e;
omegace=e*B0/me;
k=0.1/(m*v_th/e/B0);% �����Ĳ���
lambdaD=sqrt(epsilon0*T0/n0/e^2);
L=2*pi*2/k;
x0=0;
x1=L+x0;
%�������
xgrids=64;%һά�Ŀռ���
dx=L/xgrids;
xrange=x0:dx:x1; %�ռ䷶Χ
tstep=0.1/omegap;%0.5e-4*R0/v_th;%5.4143e-11;% %ʱ�䲽��
mstep=3*ceil(2*pi/omegac/tstep);
% mstep=20000;
micell=300;
zion=zeros(5,xgrids*micell); %������Ϣ 1άλ�ã�3���ٶȷ�����1��Ȩ��
vgrids=201;
vrange=[-7*v_th,7*v_th,vgrids];
dv=vrange(2)-vrange(1);
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
%% ���ӳ�ʼ��
%rng('default'); %�����������Ϊdefault��ÿ�����н����ͬ��ע�͵�����ʱ��Ϊ����
xp=linspace(x0,x1,mi)';%rand(1,mi)*(x1-x0)+x0; %����λ�ã���������ֲ�
%xp=0.008*L*cos(k*xp)+xp; %��ʼ�Ŷ�
vpx=v_th*randn(mi,1);%normrnd(0.0,v_th,[mi,1]);%�����ٶȣ����������ֲ�,ע����GTC�Ĳ��
vpy=v_th*randn(mi,1);%normrnd(0.0,v_th,[mi,1]);
vpz=v_th*randn(mi,1);%normrnd(0.0,v_th,[mi,1]);
pw=zeros(mi,1)+1;
pw=pw+0.1*sin(k*xp);
rhos=sqrt(T0*m)/e/B0;
A=zeros(xgrids,xgrids);
for i=2:xgrids-1
    A(i,i)=2/dx^2;
    A(i,i-1)=-1/dx^2;
    A(i,i+1)=-1/dx^2;
end
A(1,1)=2/dx^2;
A(1,2)=-1/dx^2;
A(1,xgrids)=-1/dx^2;
Mlap=A;
B=A;
A(end,:)=1;
A=A*n0*(m/B0^2);
A=(eye(xgrids)+rhos^2*B)*epsilon0*B+A;
B(end,1)=-1/dx^2;
B(end,end)=2/dx^2;
B(end,end-1)=-1/dx^2;
C=B;
B=B*lambdaD^2;
% for i=1:xgrids
%     B(i,i)=B(i,i)+1;
% end
B(end,:)=1;
Diff=zeros(xgrids);
for i=2:xgrids-1
    Diff(i,i+1)=1/2/dx;
    Diff(i,i-1)=-1/2/dx;
end
Diff(1,2)=1/2/dx;
Diff(1,end)=-1/2/dx;
Diff(end,1)=1/2/dx;
Diff(end,end-1)=-1/2/dx;
lambdaDe=lambdaD;
C=C/(lambdaDe^2+lambdaD^2);
for i=1:xgrids
    C(i,i)=C(i,i)+1/lambdaD^2;
end
C(end,:)=1;


%%
figure;
dtime=tstep;
for istep=1:mstep
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
%     xp=mod(xp-x0,(x1-x0))+x0;%�������꣬ʹ���������Ա߽�����
    
    g0=floor(xp/dx-.5)+1;%floor((xp-x0)/dx-0.5)+1;
    g=[g0;g0+1];
    out=(g<1);g(out)=g(out)+xgrids;
    out=(g>xgrids);g(out)=g(out)-xgrids;
    h1=abs(xp/dx-g0+.5);%mod((xp-x0)/dx+0.5,1); 
    h=[1-h1;h1];
    mat=sparse(parray,g,h,mi,xgrids);
    pw=pw*(n0/micell);
    rho=full(e/epsilon0*(pw'*mat))-rho0;
    pw=pw/(n0/micell);
%     rho=0.05*rho0*cos(k*[dx/2:dx:L-dx/2]);

%     poissonE(:)=0.0;
%     for xi=1:xgrids-1
%         poissonE(xi+1)=poissonE(xi)+0.5*(rho(xi)+rho(xi+1))*dx;
%     end
%     poissonE(1)=poissonE(xgrids)+0.5*(rho(xgrids)+rho(1))*dx;
%     poissonE=poissonE-mean(poissonE);
% %     poissonE=poissonE/(1+omegape^2/omegace^2);
%     phi=zeros(xgrids,1);
%     for xi=1:xgrids-1
%         phi(xi+1)=phi(xi)-0.5*(poissonE(xi)+poissonE(xi+1))*dx;
%     end
%     phi(1)=phi(xgrids)-0.5*(poissonE(1)+poissonE(xgrids))*dx;
%     phi=phi-mean(phi);
%     poissonE=poissonE+extE';

    phi=B\(lambdaD^2*rho');
    poissonE=-Diff*phi;
%     for xi=2:xgrids-1
%         poissonE(xi)=(-phi(xi+1)+phi(xi-1))/2/dx;
%     end
%     poissonE(1)=(-phi(2)+phi(end))/2/dx;
%     poissonE(end)=(-phi(1)+phi(end-1))/2/dx;
%     poissonE=poissonE-mean(poissonE);


    vpx=vpx+mat*poissonE*e/m*dtime/2;
    gammab=e*B0/m*dtime/(1+e^2*B0^2/m^2/4*dtime^2);
    vx1=vpx+(vpy-vpx*e*B0/m*dtime/2)*gammab;
    vpy=vpy-(vpx+vpy*e*B0/m*dtime/2)*gammab;
    vpx=vx1;
    vpx=vx1+mat*poissonE*e/m*dtime/2;
    

%     vpx=vpx+e/m*(mat*poissonE-vpz*B0)*dtime;
%     vpy=vpy-e/m*(vpx*B0)*dtime;
%     vpz=vpz+e/m*(vpx*B0)*dtime;
    %pw=pw+e/T0*mat*poissonE.*vpx*dtime;
    
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
    pw=pw*(n0/micell);
    g_rho=(full(e/epsilon0*(pw'*mat))-rho0)*epsilon0*0.5...
          +0.25*(full(e/epsilon0*(pw'*mat_left))-rho0)*epsilon0...
          +0.25*(full(e/epsilon0*(pw'*mat_right))-rho0)*epsilon0;
    g_rho=g_rho-mean(g_rho);
    pw=pw/(n0/micell);
    %g_rho=0.05*rho0*epsilon0*cos(k*[dx/2:dx:L-dx/2]);
    g_rho1=g_rho';
%     g_rho1=(eye(xgrids)-rhos^2*Mlap)*g_rho1;
    g_rho1(end)=0;
%     g_rho1=g_rho1;
    g_phi=A\g_rho1;
%     g_phi=g_phi+g_rho'/(n0*e^2/T0);
    %g_phi=g_phi/(1+omegape^2/omegace^2);
    
    %diagnosis
%     if istep>mstep
    Etime(:,istep)=phi;
    gyro_phi(:,istep)=g_phi;
%     Etime(:,istep)=poissonE;
    ntime(:,istep)=rho;
%     end

    subplot(3,2,1);
    plot(dx/2:dx:L-dx/2,poissonE);
    title(['t=',num2str(istep)]);
    %ylim([-1e6,1e6]);
    subplot(3,2,3);plot(dx/2:dx:L-dx/2,rho*epsilon0/e/n0);
    subplot(3,2,4);plot(dx/2:dx:L-dx/2,g_rho/e/n0);
    subplot(3,2,5);plot(dx/2:dx:L-dx/2,phi);
    subplot(3,2,6);plot(dx/2:dx:L-dx/2,g_phi);
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
mean_n=mean(ntime,2);
subplot(211);
plot(mean_E);
hold on;
plot(g_phi(:,end),'r--'); hold off;
subplot(212);
plot(Etime(xgrids/2,:));

% subplot(212);plot(mean_n*epsilon0/e/n0);hold on; plot(g_rho/e/n0,'r--');hold off;

