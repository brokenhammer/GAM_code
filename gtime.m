function dgdt = gtime(g,phi,omega_di,q,R0,pgpth,kr,vperp,omega0,vpara)
%DGDT calc time derivative of g using phi
%   此处显示详细说明
pgpth(2:end-1,:,:) = (g(3:end,:,:) - g(1:end-2,:,:))/2/dtheta;
pgpth(1,:,:) = (g(2,:,:) - g(end-1,:,:))/2/dtheta;
pgpth(end,:,:) = pgpth(1,:,:);
krrho = kr * vperp / omega0;

dgdt = -vpara/q/R0.*pgpth - 1i*omega_di.*g - 1i * omega_di.*besselj(0,krrho).*Fm.*phi*e/Ti;
end

