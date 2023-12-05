%% Compensateurs inclinaison de la plaque (Cas 2)

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

% Conception de l'AvPh à partir du diagramme de bode
[Gm_k,Pm_k,Wp_k,Wg_k] = margin(FT_K_desire);
avphreq = (pm-Pm_k)/2;
avphreq_5 = avphreq + 5;
alpha = (1-sind(avphreq_5))/(1+sind(avphreq_5));
T = 1/(wg*sqrt(alpha));
num = conv([1 1/T],[1 1/T]);
den = conv([1 1/(alpha*T)],[1 1/(alpha*T)]);
G_AvPh_0gain = tf(num,den);
ka = (K_desire/alpha)- 13395;
G_AvPh = ka*G_AvPh_0gain;
FT_AvPh = FTBO*G_AvPh;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Conception d'un PI %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[num,den] = tfdata(FT_AvPh,'v');
K_pos_actuel = num(end)/den(end);
beta = (kpos_desire/K_pos_actuel)-1;
T_RePh = 4/wg;
z = -wg/10;
G_PI = tf([1 -z],[1 0]);
FT_DAvPh_PI = FT_AvPh*G_PI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Affichage des graphiques %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphique de l'échelon
t = (0:0.0001:0.5);
u = ones(size(t));
y = lsim(feedback(FT_DAvPh_PI,1),u,t);
figure(1)
plot(t,y)

% Lieu de Bode
figure(2)
margin(FT_DAvPh_PI)

% Lieu de Nyquist
figure(3)
nyquist(FT_DAvPh_PI);

% Lieu des racines 
figure(4)
rlocus(FT_DAvPh_PI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


