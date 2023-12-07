%% Compensateurs inclinaison de la plaque (Cas 2)

%%%%%%%%%%%%%%%%%%% Commande de base %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Spécifications à atteindre %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pm_hpc2 = 25;
wg_hpc2 = 185;
consigne_echelon1_hpc2 = 0.01;
erreur_echelon1_hpc2 = -0.000425;
erreur_echelon2_hpc2 = 0;
FTBO_hpc2 = tf(29.36,[1 31.3 -1216 -3.805e04]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Conception d'un Double AvPh %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zeta = 0.5*sqrt(tand(pm)*sind(pm))
% wn = (wg*tand(pm))/(2*zeta)
kpos_desire_hpc2 = 1/erreur_echelon2_hpc2;
[mag_hpc2,pha_hpc2]=bode(FTBO_hpc2,wg_hpc2);
K_desire_hpc2 = 1/mag_hpc2;

% Fonction de transfert
FT_K_desire_hpc2 = FTBO_hpc2*K_desire_hpc2;

% Conception de l'AvPh à partir du diagramme de bode
[Gm_k_hpc2,Pm_k_hpc2,Wp_k_hpc2,Wg_k_hpc2] = margin(FT_K_desire_hpc2);
avphreq_hpc2 = (pm_hpc2-Pm_k_hpc2)/2;
avphreq_5_hpc2 = avphreq_hpc2 + 5;
alpha_hpc2 = (1-sind(avphreq_5_hpc2))/(1+sind(avphreq_5_hpc2));
T_hpc2 = 1/(wg_hpc2*sqrt(alpha_hpc2));
numDAvPh_hpc2 = conv([1 1/T_hpc2],[1 1/T_hpc2]);
denDAvPh_hpc2 = conv([1 1/(alpha_hpc2*T_hpc2)],[1 1/(alpha_hpc2*T_hpc2)]);
G_1_hpc2 = tf(numDAvPh_hpc2,denDAvPh_hpc2);
ka_hpc2 = (K_desire_hpc2/alpha_hpc2)- 13395;
G_DAvPh_hpc2 = ka_hpc2*G_1_hpc2;
FT_DAvPh_hpc2 = FTBO_hpc2*G_DAvPh_hpc2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Conception d'un PI %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[numDAvPh_hpc2,denDAvPh_hpc2] = tfdata(FT_DAvPh_hpc2,'v');
K_pos_actuel_hpc2 = numDAvPh_hpc2(end)/denDAvPh_hpc2(end);
beta_hpc2 = (kpos_desire_hpc2/K_pos_actuel_hpc2)-1;
T_RePh_hpc2 = 4/wg_hpc2;
z_hpc2 = -wg_hpc2/10;
G_PI_hpc2 = tf([1 -z_hpc2],[1 0]);
G_hpc2 = G_PI_hpc2*G_DAvPh_hpc2;
FT_DAvPh_PI_hpc2 = FT_DAvPh_hpc2*G_PI_hpc2;
[num_hpc2,den_hpc2] = tfdata(G_hpc2,'v'); % Pour SIMULINK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Affichage des graphiques %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphique de l'échelon
t = (0:0.0001:0.5);
u = ones(size(t));
y = lsim(feedback(FT_DAvPh_PI_hpc2,1),u,t);
figure(1)
plot(t,y)

% Lieu de Bode
figure(2)
margin(FT_DAvPh_PI_hpc2)

% Lieu de Nyquist
figure(3)
nyquist(FT_DAvPh_PI_hpc2);

% Lieu des racines 
figure(4)
rlocus(FT_DAvPh_PI_hpc2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


