%%
% @file: Capteur.m
% @Date: 2023/09/29
%
% @Name: - Shawn Miller (mils2203)
% @Brief: equation de tension par rapport a la distance
%%


close all
clc

% d = K*exp(v/vs)
% log(d) = log(K) + v/vs
% D = log(d)
% V = v
% m = -1/vs
% b = log(K)



distance(end) = distance(end-1);
N = length(distance)-1;

phiM = [N sum(voltage);sum(voltage) sum(voltage.^2)];
D = [sum(distance);sum(distance.*voltage)];

coeffs = inv(phiM)*D;
vs = coeffs(2)
K = coeffs(1)

D_moy = mean(log(distance));
D_v = log(K) + voltage./vs;
d_app = K*exp((voltage./vs));
%d_app = K + vs.*voltage;
%D_v = d_app;

E_RMS = sqrt(mean((D_v-log(distance)).^2))
E_R2  = sum((D_v-D_moy).^2)/sum((log(distance)-D_moy).^2)

figure()
plot(voltage,distance)
hold on
plot(voltage,d_app, '-')
grid on
ylabel("Distance plaque (m)")
xlabel("Tension (V)")
title("Distance D selon Tension A")




