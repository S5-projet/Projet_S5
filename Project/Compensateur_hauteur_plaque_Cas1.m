%% Compensateurs inclinaison de la plaque (Cas 1)

%%%%%%%%%%%%%%%%%%% Commande de base %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear all
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Spécifications à atteindre %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pm_hpc1 = 25;
wg_hpc1 = 185;
consigne_echelon1_hpc1 = 0.01;
erreur_echelon1_hpc1 = -0.000425;
erreur_echelon2_hpc1 = 0;
FTBO_hpc1 = tf(29.36,[1 31.3 -1216 -3.805e04]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Conception d'un Double AvPh %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta = 0.5*sqrt(tand(pm_hpc1)*sind(pm_hpc1));
wn = (wg_hpc1*tand(pm_hpc1))/(2*zeta);
bw = wn*sqrt((1-2*zeta^2)+sqrt(4*zeta^4-4*zeta^2+2));
kpos_desire_hpc1 = 1/erreur_echelon1_hpc1;
[mag_hpc1,pha_hpc1]=bode(FTBO_hpc1,wg_hpc1);
K_desire_hpc1 = 1/mag_hpc1;

% Fonction de transfert
FT_K_desire_hpc1 = FTBO_hpc1*K_desire_hpc1;

[Gm_k_hpc1,Pm_k_hpc1,Wp_k_hpc1,Wg_k_hpc1] = margin(FT_K_desire_hpc1);
avphreq_hpc1 = (pm_hpc1-Pm_k_hpc1)/2;
avphreq_5_hpc1 = avphreq_hpc1 + 7.5;
alpha_hpc1 = (1-sind(avphreq_5_hpc1))/(1+sind(avphreq_5_hpc1));
T_hpc1 = 1/(wg_hpc1*sqrt(alpha_hpc1));
numDAvPh_hpc1 = conv([1 1/T_hpc1],[1 1/T_hpc1]);
denDAvPh_hpc1 = conv([1 1/(alpha_hpc1*T_hpc1)],[1 1/(alpha_hpc1*T_hpc1)]);
G1_hpc1 = tf(numDAvPh_hpc1,denDAvPh_hpc1);
ka_hpc1 = (K_desire_hpc1/alpha_hpc1)-95350;
G_DAvPh_hpc1 = ka_hpc1*G1_hpc1;
FT_DAvPh_hpc1 = FTBO_hpc1*G_DAvPh_hpc1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Conception d'un RePh %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[numDAvPh_hpc1,denDAvPh_hpc1] = tfdata(FT_DAvPh_hpc1,'v');
K_pos_actuel_hpc1 = numDAvPh_hpc1(end)/denDAvPh_hpc1(end);
beta_hpc1 = (kpos_desire_hpc1/K_pos_actuel_hpc1)-1;
T_RePh_hpc1 = 4/wg_hpc1;
z_hpc1 = -1/T_RePh_hpc1;
p_hpc1 = -1/(beta_hpc1*T_RePh_hpc1);
G_RePh_hpc1 = tf([1 -z_hpc1],[1 -p_hpc1]);
FT_DAvPh_RePh_hpc1 = FT_DAvPh_hpc1*G_RePh_hpc1;
G_hpc1 = G_DAvPh_hpc1*G_RePh_hpc1;
[num_hpc1,den_hpc1] = tfdata(G_hpc1,'v') ; % Pour SIMULINK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Affichage des graphiques %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphique de l'échelon
t_hpc1 = (0:0.0001:0.5);
u_hpc1 = ones(size(t_hpc1));
y_hpc1 = lsim(feedback(FT_DAvPh_RePh_hpc1,1),u_hpc1,t_hpc1);

figure(1)
plot(t_hpc1,y_hpc1*-consigne_echelon1_hpc1)
title('Réponse à un échelon de 0.01')
xlabel('Temps (s)')
ylabel('Hauteur de la plaque')
grid minor
yline(-0.0104,'r')


% Lieu de Bode
figure(2)
margin(FTBO_hpc1)
hold on
margin(FT_DAvPh_RePh_hpc1);
legend('FT','FT compensé')

% Lieu de Nyquist
figure(3)
nyquist(FT_DAvPh_RePh_hpc1);

% Lieu des racines 
figure(4)
rlocus(FTBO_hpc1,'b');
hold on
rlocus(FT_DAvPh_RePh_hpc1,'r');
legend('FT','FT compensé')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testdiscret(G_hpc1)

