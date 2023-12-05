%% Compensateurs inclinaison de la plaque (Cas 1)

%%%%%%%%%%%%%%%%%%% Commande de base %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Spécifications à atteindre %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pm = 25;
wg = 185;
consigne_echelon1 = 0.01;
erreur_echelon1 = -0.000425;
erreur_echelon2 = 0;
FTBO = tf(29.36,[1 31.3 -1216 -3.805e04]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Conception d'un Double AvPh %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zeta = 0.5*sqrt(tand(pm)*sind(pm))
% wn = (wg*tand(pm))/(2*zeta)
kpos_desire = 1/erreur_echelon2;
[mag,pha]=bode(FTBO,wg);
K_desire = 1/mag;

% Fonction de transfert
FT_K_desire = FTBO*K_desire;

[Gm_k,Pm_k,Wp_k,Wg_k] = margin(FT_K_desire);
avphreq = (pm-Pm_k)/2;
avphreq_5 = avphreq + 7.5;
alpha = (1-sind(avphreq_5))/(1+sind(avphreq_5));
T = 1/(wg*sqrt(alpha));
num = conv([1 1/T],[1 1/T]);
den = conv([1 1/(alpha*T)],[1 1/(alpha*T)]);
G_AvPh_0gain = tf(num,den);
ka = (K_desire/alpha)-95350;
G_AvPh = ka*G_AvPh_0gain;
FT_AvPh = FTBO*G_AvPh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Conception d'un RePh %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num,den] = tfdata(FT_AvPh,'v');
K_pos_actuel = num(end)/den(end);
beta = (kpos_desire/K_pos_actuel)-1;
T_RePh = 4/wg;
z = -1/T_RePh;
p = -1/(beta*T_RePh);
G_RePh = tf([1 -z],[1 -p]);
FT_DAvPh_RePh = FT_AvPh*G_RePh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Affichage des graphiques %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphique de l'échelon
t = (0:0.0001:0.5);
u = ones(size(t));
y = lsim(feedback(FT_DAvPh_RePh,1),u*consigne_echelon1,t);
figure(1)
plot(t,y)

% Lieu de Bode
figure(2)
margin(FT_DAvPh_RePh)

% Lieu de Nyquist
figure(3)
nyquist(FT_DAvPh_RePh);

% Lieu des racines 
figure(4)
rlocus(FT_DAvPh_RePh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


