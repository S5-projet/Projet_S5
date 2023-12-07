%% Compensateurs inclinaison de la plaque

%%%%%%%%%%%%%%%%%%% Commande de base %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear all
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Spécifications à atteindre %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mp_ip = 5;
ts_ip = 0.03;
tp_ip = 0.025;
tr_ip = 0.02;
Erreur_ip = 0;
FTBO_ip = tf(9633,[1 31.3 -1807 -5.658e04]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Conception d'un Double AvPh %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_ip = atand(-pi/log(mp_ip/100));
zeta_ip = cosd(phi_ip);
wn1_ip = 4/(zeta_ip*ts_ip);
wn2_ip = pi/(tp_ip*sqrt(1-zeta_ip^2));
wn3_ip = (1+1.1*zeta_ip+1.4*zeta_ip^2)/tr_ip;
wn_ip = wn1_ip;
bw = wn_ip*sqrt((1-2*zeta_ip^2)+sqrt(4*zeta_ip^4-4*zeta_ip^2+2));
pole_p_ip = (-zeta_ip*wn_ip)+(i*(wn_ip*(sqrt(1-zeta_ip^2))));
pole_n_ip = (-zeta_ip*wn_ip)-(i*(wn_ip*(sqrt(1-zeta_ip^2))));

phase_pole_des_ip = rad2deg(angle(evalfr(FTBO_ip,pole_p_ip))) - 360;
delta_phi_ip = ((-180 - phase_pole_des_ip)/2)+0.05;  %%+8.45;
alpha_ip = 180 - phi_ip;

phi_z_ip = (alpha_ip+delta_phi_ip)/2;
phi_p_ip = (alpha_ip-delta_phi_ip)/2;

z_ip = real(pole_p_ip) - (imag(pole_p_ip)/tand(phi_z_ip));
p_ip = real(pole_p_ip) - (imag(pole_p_ip)/tand(phi_p_ip));

G1_ip = (tf([1 -z_ip],[1 -p_ip]))^2;
FT1_ip = FTBO_ip * G1_ip;
ka_ip = (1/ abs(evalfr(FT1_ip,pole_p_ip)))+12500;
G_DAvPh_ip = ka_ip*G1_ip;
FT_DAvPh_ip = FTBO_ip*G_DAvPh_ip;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Conception d'un PI %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_pi_ip = real(pole_p_ip)/1.75;
ki_ip = 1;
G_PI_ip = tf([1 -z_pi_ip],[1 0]);
FT_DAvPh_PI_ip = FT_DAvPh_ip * G_PI_ip;
G_ip = G_DAvPh_ip*G_PI_ip;
[num_ip,den_ip] = tfdata(G_ip,'v'); % Pour SIMULINK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Affichage des graphiques %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphique de la réponse à l'échelon
t_ip = (0:0.0001:0.5);
u_ip = ones(size(t_ip));
y_ip = lsim(feedback(FT_DAvPh_PI_ip,1),u_ip,t_ip);

figure(1)
plot(t_ip,y_ip)
hold on
title('Réponse à un échelon unitaire')
xlabel('Temps (s)')
ylabel('Angle d''inclinaison')
grid minor
yline(1.02,'r')
yline(0.98,'r')

% Lieu de Bode
figure(2)
margin(FTBO_ip)
hold on
margin(FT_DAvPh_PI_ip);
legend('FT','FT compensé')


% Lieu de Nyquist
figure(3)
nyquist(FT_DAvPh_PI_ip);

% Lieu des racines 
figure(4)
rlocus(FTBO_ip,'b')
hold on
rlocus(FT_DAvPh_PI_ip,'r');
legend('FT','FT compensé')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testdiscret(G_ip)


