%% 
% @file: run_bancessai.m
% @Date: 2023/09/29
%
% @Name: - Shawn Miller (mils2203)
%        -           
%        -   
%        -   
%        -   
%        -   
%        -   
%
% @Brief: 
%   
%
%%

close all
clear all
clc

% Position à l'équilibre de la sphère (pour tests statiques)
sig = 1.0;         % Présence (1) ou non (0) de la sphère
xSeq = 0.000;      % Position x de la sphère à l'équilibre en metres
ySeq = 0.000;      % Position y de la sphère à l'équilibre en metres

%Point d'opération choisi pour la plaque
Axeq = 0;               %en degres
Ayeq = 0;               %en degres
Pzeq = .015;            %en metres
Value_initial = 1;

%Exemple de trajectoire
t_des     = [0:1:8]'*5;
x_des     = [t_des, [0 0 0.5 1  0 -1 0 1 0]'*0.05];
y_des     = [t_des, [0 0 0 0 -1  0 1 0 0]'*0.05];
z_des     = [t_des, [1 1 1 1  1  1 1 1 1]'*.015];
tfin = 50;

%initialisation
bancessai_ini  %faites tous vos calculs de modele ici
bancEssaiConstantes

%Calcul des compensateurs
%iniCTL_ver4    %Calculez vos compensateurs ici

%simulation
open_system('DYNctl_ver4_etud_obfusc')
set_param('DYNctl_ver4_etud_obfusc','AlgebraicLoopSolver','LineSearch')
sim('DYNctl_ver4_etud_obfusc')

%affichage

%trajectoires
figure(1)
subplot(1,2,1)
hold on;
axis([-0.06 0.06 -0.06 0.06]);
plot(ynonlineaire(:,7),ynonlineaire(:,8));
title('Trajectoire attendue')
subplot(1,2,2)
axis([-0.06 0.06 -0.06 0.06]);
plot(ynonlineaire1(:,8),ynonlineaire1(:,7));
title('Trajectoire effectué')
hold off;

figure(4)
hold on;
plot(ylineaire(:,4),ylineaire(:,5));
title('Trajectoire du système linéaire')

figure()
 hold on;
 plot(yDecouple(:,4),yDecouple(:,5));
 title('Trajectoire du système découplé')


