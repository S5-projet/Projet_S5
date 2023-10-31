%% Commande de base

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ajout de précision lors de l'affichage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chargement des données des capteurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Fe_attraction.mat');
load('Fs.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Coefficients de Fs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Fs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs_1 = ones(size(z_pos));
Fs_X = [Fs_1, z_pos.^1, z_pos.^2, z_pos.^3, z_pos.^4, z_pos.^5, z_pos.^6]; 
Fs_Y = Fs;
Fs_coeff = pinv(Fs_X)*(-1./Fs_Y);
Fs_coeff2 = -1./(Fs_coeff(1) + Fs_coeff(2).*z_pos + Fs_coeff(3).*z_pos.^2 + Fs_coeff(4).*z_pos.^3 + Fs_coeff(5).*z_pos.^4 + Fs_coeff(6).*z_pos.^5 + Fs_coeff(7).*z_pos.^6);
%disp('Coefficients de Fs:');
%disp(Fs_coeff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Graphique de comparaison du Fs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1);
% plot(z_pos,Fs_coeff2);
% title('Lissage du Fs');
% hold on;
% plot(z_pos,Fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Coefficients de Fek

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Fe_m1A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fe_m1A_1 = ones(size(z_m1A));
Fe_m1A_X = [Fe_m1A_1, z_m1A.^1, z_m1A.^2, z_m1A.^3]; 
Fe_m1A_Y = Fe_m1A;

be = 13.029359254409743;
ik_1 = -1;
num_m1A = (ik_1^2+be*abs(ik_1))*sign(ik_1);
Fe_m1A_coeff = pinv(Fe_m1A_X)*(num_m1A./Fe_m1A_Y);
Fe_m1A_coeff2 = num_m1A./(Fe_m1A_coeff(1) + Fe_m1A_coeff(2).*z_m1A + Fe_m1A_coeff(3).*z_m1A.^2 + Fe_m1A_coeff(4).*z_m1A.^3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Graphique de comparaison du Fe_m1A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2);
% plot(z_m1A,Fe_m1A_coeff2);
% title('Lissage du Fe-m1A');
% hold on;
% plot(z_m1A,Fe_m1A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Fe_m2A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fe_m2A_1 = ones(size(z_m2A));
Fe_m2A_X = [Fe_m2A_1, z_m2A.^1, z_m2A.^2, z_m2A.^3]; 
Fe_m2A_Y = Fe_m2A;

be = 13.029359254409743;
ik_2 = -2;
num_m2A = (ik_2^2+be*abs(ik_2))*sign(ik_2);
Fe_m2A_coeff = pinv(Fe_m2A_X)*(num_m2A./Fe_m2A_Y);
Fe_m2A_coeff2 = num_m2A./(Fe_m2A_coeff(1) + Fe_m2A_coeff(2).*z_m2A + Fe_m2A_coeff(3).*z_m2A.^2 + Fe_m2A_coeff(4).*z_m2A.^3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Graphique de comparaison du Fe_m2A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3);
% plot(z_m2A,Fe_m2A_coeff2);
% title('Lissage du Fe-m2A');
% hold on;
% plot(z_m2A,Fe_m2A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Moyenne des coefficients du Fe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fe_coeff = (Fe_m1A_coeff + Fe_m2A_coeff)/2;
%disp('Coefficients de Fe:');
%disp(Fe_coeff);
num_fe = (num_m1A + num_m2A)/2;
z_fe = (z_m1A(1:end-5) + z_m2A)/2;
Fe = (Fe_m1A(1:end-5)+Fe_m2A)/2;
Fe_coeff2 = num_fe./(Fe_coeff(1) + Fe_coeff(2).*z_fe + Fe_coeff(3).*z_fe.^2 + Fe_coeff(4).*z_fe.^3);

%% Calcul de l'erreur RMS et de la corrélation R^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Erreur RMS et Corrélation R^2 de Fs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erreur RMS
N = 7;
E_Fs = sum((Fs_coeff2-Fs).^2);
RMS_Fs = sqrt((1/N)*E_Fs);
%disp('Erreur RMS de Fs:');
%disp(RMS_Fs);

% Coefficient de corrélation R^2
ybarre = (1/N)*(sum(Fs));
r_Fs = sum((Fs_coeff2-ybarre).^2)./sum((Fs - ybarre).^2);
%disp('Coefficient de corrélation R^2 de Fs:');
%disp(r_Fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Erreur RMS et Corrélation R^2 de Fe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erreur RMS
N = 4;
E_Fe = sum((Fe_coeff2-Fe).^2);
RMS_Fe = sqrt((1/N)*E_Fe);
%disp('Erreur RMS de Fe:');
%disp(RMS_Fe);

% Coefficient de corrélation R^2
ybarre = (1/N)*(sum(Fe));
r_Fe = sum((Fe_coeff2-ybarre).^2)./sum((Fe - ybarre).^2);
%disp('Coefficient de corrélation R^2 de Fe:');
%disp(r_Fe);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
