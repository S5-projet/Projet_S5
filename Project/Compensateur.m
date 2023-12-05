%% Compensateur position de la sphère

clc
clear all
close all

FT_sx_p = tf([-7.007],[1 0 0]);
FT_sy_p = tf([7.007],[1 0 0]);

[num_pos_sph,den_pos_sph] = tfdata(FT_sx_p,'v');

ts_posi_sphere = 2;
zeta_position_sphere = 0.9;
wn_position_sphere = 4/(zeta_position_sphere*ts_posi_sphere);

pole1 = -zeta_position_sphere*wn_position_sphere + i*wn_position_sphere*sqrt(1-(zeta_position_sphere));
pole2 = -zeta_position_sphere*wn_position_sphere - i*wn_position_sphere*sqrt(1-(zeta_position_sphere));

% newGs_position_sphere = 7.007/(pole1^2);
% zPI_ = real(pole1)/10;
% Kp = 1/(abs(((pole1-zPI)/pole1)*newGs_position_sphere));
% CompensateurPI = Kp*tf([pole1 -zPI],[pole1]);

A = num_pos_sph(3);%7.007
Kv = (2*zeta_position_sphere*wn_position_sphere)/A;
Kp = (wn_position_sphere^2)/A;

FT_sx_p_asservie = tf([Kp*A],[1 (A*Kv) (A*Kp)]);

figure(1)
rlocus(FT_sx_p);

figure(2)
margin(FT_sx_p);

figure(3)
rlocus(FT_sx_p_asservie);

figure(4)
margin(FT_sx_p_asservie);