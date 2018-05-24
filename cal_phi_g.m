function phi = cal_phi_g(g, vperp, omega0, dvpe, dvpa, kr)
%CAL_PHI_G 从delta_g 积分计算 delta phi， 0th Maxwellian distribution.
%   此处显示详细说明

krrho = kr * vperp / omega0;
deltan = 2*pi*dvpe*dvpa * vperp.*g.*besselj(0,krrho);
deltan = sum(deltan,2);
deltan = squeeze(sum(deltan,3));
    
phi = deltan ./ (e * n0 / Ti * (1-gamma0));

