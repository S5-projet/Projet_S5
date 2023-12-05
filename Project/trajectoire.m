clc 
close all
clear all
% x=[6 5 4 3 2 1 0]';
% y=[8 2 3 3 7 1 0]';

x=[0 1 2 3 4 5 6]';
y=[0 1 7 3 3 2 8]';

v= 10;
ts= 0.04;

[Pi,Ltr, E, Vr, Traj, tt] = traject(x,y,v,ts);

figure
hold on
plot(Traj(:, 1), Traj(:, 2), ':xk');
scatter(x, y, 'b', 'linewidth', 2);
hold off

[N M] = size(x);
[L K] = size(Traj);
dist = zeros(L, 1);
index = zeros(N, 1);
for i = 1:N
    for j = 1:L
        dist(j) = sqrt((y(i) - Traj(j, 2))^2 + (x(i) - Traj(j, 1))^2);
    end
    %disp(dist)
    [bogus index(i)] = min(dist);
end

for i = 1:N
    Traj_E(i, :) = [Traj(index(i), 1) Traj(index(i), 2)];
end
R_carre = sum(Traj_E(:,2)-(mean(y))^2)/sum(y-(mean(y))^2);
err_RMS = sqrt(mean((Traj_E(:,2)-y).*(Traj_E(:,2)-y)));
disp("R2 : " + R_carre)
disp("RMS : " + err_RMS)

function [Pi,Ltr, E, Vr, Traj, tt]= traject(x,y,v,ts)

    % polynome d'interpolation
    
    N = length(x);
    
    for i=1:N
        roger(:,i)= x.^(i-1);
    end 
    
    Pi= pinv(roger)*y;
    Pi=flip(Pi);
    
    
    %Calcul longueur trajectoire
    
    M=101;
    dx = linspace(x(1), x(end), M);

    fy=polyval(Pi,dx);

    figure
    plot(x,y)
    hold on
    plot(dx, fy)

    
    for i=1:length(Pi)-1
        Pi_d(i)= (N-i)*Pi(i);
    end 
    
    yd= polyval(Pi_d,dx);
    
    
    for j=1:length(Pi_d)-1
        ddpi(j)= (N-j)*Pi_d(j);
    end 
    
    ydd= polyval(ddpi,dx);
    
    g= sqrt(1+(yd).^2);
    
    h= abs((dx(2)-dx(1)));
    
    for i=2:M
        Ltr(i) = (h/2)*(g(1) + g(i) + 2*sum(g(2:1:i-1)));
    end
    
    %calcul de l'erreur 
    dg = (yd.*ydd./g);
    
    E = (h^2/12)*(((dg(end)-dg(end-1))/h)-((dg(2)-dg(1))/h));

    
    % calcul de O et vitesse réel 
    
    d = v*ts;
    O = round(Ltr(end)/(d));
    Vr = Ltr(end)/(ts*(O));
    tt = Ltr(end)/Vr;
    dl= Ltr(end)/O;
    % Newton-Raphson
    
    %disp(O)
    fn=1000;
    an=min(x);
    bn = an + dl;
    X(1)=an;
    it=0;
    tol=1e-08;
    for i=1:O
        while abs(fn) > tol && it < 1000
            xn = linspace(an, bn, 15);
            yd= polyval(Pi_d,xn);
            g= sqrt(1+(yd).^2);
            dn= g(end);
            h1=(xn(2)-xn(1));
            fn = ((g(1) + g(15) + 4*sum(g(2:2:15-1)) + 2*sum(g(3:2:15-1)))*h1/3)-dl;
            bn = bn - ((fn)/dn);
            it=it+1;
            %disp("an : " + an + ", bn : " + bn)
        end 
        it=0;
        X(i)=bn;
        an = bn;
        fn=0.01;
        
    end
    Y = polyval(Pi,X);
    %disp("O" + O)
    Traj = [X' Y'];
   
end 
 


