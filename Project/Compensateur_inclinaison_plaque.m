%% Compensateurs inclinaison de la plaque

%%%%%%%%%%%%%%%%%%% Commande de base %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Spécifications à atteindre %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mp = 5;
ts = 0.03;
tp = 0.025;
tr = 0.02;
Erreur = 0;
FTBO = tf(9633,[1 31.3 -1807 -5.658e04]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Conception d'un Double AvPh %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = atand(-pi/log(mp/100));
zeta = cosd(phi);
wn1 = 4/(zeta*ts);
wn2 = pi/(tp*sqrt(1-zeta^2));
wn3 = (1+1.1*zeta+1.4*zeta^2)/tr;
wn = wn1;
pole_p = (-zeta*wn)+(i*(wn*(sqrt(1-zeta^2))));
pole_n = (-zeta*wn)-(i*(wn*(sqrt(1-zeta^2))));

phase_pole_des = rad2deg(angle(evalfr(FTBO,pole_p))) - 360;
delta_phi = ((-180 - phase_pole_des)/2)+8.45;
alpha = 180 - phi;

phi_z = (alpha+delta_phi)/2;
phi_p = (alpha-delta_phi)/2;

z = real(pole_p) - (imag(pole_p)/tand(phi_z));
p = real(pole_p) - (imag(pole_p)/tand(phi_p));

G_AvPh_0gain = (tf([1 -z],[1 -p]))^2;
FTBO_AvPH_0gain = FTBO * G_AvPh_0gain;
ka = (1/ abs(evalfr(FTBO_AvPH_0gain,pole_p)))+12500;
G_AvPh_gain = ka*G_AvPh_0gain;
FT_AvPh_gain = FTBO*G_AvPh_gain;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Conception d'un PI %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_pi = real(pole_p)/1.75;
ki = 1;
G_PI = tf([1 -z_pi],[1 0]);
FT_DAvPh_PI = FT_AvPh_gain*G_PI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Affichage des graphiques %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphique de la réponse à l'échelon
t = (0:0.0001:0.5);
u = ones(size(t));
y = lsim(feedback(FT_DAvPh_PI,1),u,t);
figure(1)
plot(t,y)

% Lieu de Bode
figure(2)
margin(FT_DAvPh_PI);

% Lieu de Nyquist
figure(3)
nyquist(FT_DAvPh_PI);

% Lieu des racines 
figure(4)
rlocus(FT_DAvPh_PI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


