%% Compensateur position de la sphère

FT_sx_p = tf([-7.007],[1 0 0]);
FT_sy_p = tf([7.007],[1 0 0]);

[num_pos_sph,den_pos_sph] = tfdata(FT_sx_p,'v');

ts_pos_sph = 2;
zeta_pos_sph = 0.9;
wn_pos_sph = 4/(zeta_pos_sph*ts_pos_sph);

pole1_pos_sph = -zeta_pos_sph*wn_pos_sph + i*wn_pos_sph*sqrt(1-(zeta_pos_sph));
pole2_pos_sph = -zeta_pos_sph*wn_pos_sph - i*wn_pos_sph*sqrt(1-(zeta_pos_sph));

A_pos_sph = num_pos_sph(3);%7.007
Kv_pos_sph = (2*zeta_pos_sph*wn_pos_sph)/A_pos_sph;
Kp_pos_sph = (wn_pos_sph^2)/A_pos_sph;

FT_sx_p_asservie = tf([Kp_pos_sph*A_pos_sph],[1 (A_pos_sph*Kv_pos_sph) (A_pos_sph*Kp_pos_sph)]);
[num_comp_pos_sph,den_comp_pos_sph] = tfdata(FT_sx_p_asservie,'v');

figure(1)
rlocus(FT_sx_p);

figure(2)
margin(FT_sx_p);

figure(3)
rlocus(FT_sx_p_asservie);

figure(4)
margin(FT_sx_p_asservie);