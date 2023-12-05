clc 
close all 
<<<<<<< HEAD
clear all
x=[0 1 2 3 4 5 6]';
y=[0 1 7 3 3 2 8]';
=======

x=[0 1 7 3 3 2 8]'
y=[0 1 2 3 4 5 6]'
>>>>>>> ee520c7bf0def907087a5587872453e16578686d

v= 0.5;
ts= 0.04;

traject(x,y,v,ts);
function [Pi,Ltr, E, Vr, Traj, tt, tab]= traject(x,y,v,ts)

    % polynome d'interpolation
    
    N = length(x);
    
    for i=1:N
        roger(:,i)= x.^(i-1);
    end 
    
    Pi= pinv(roger)*y;
    Pi=flip(Pi);
    
    
    %Calcul longueur trajectoire
    
    M=15;
    dx = linspace(x(1), x(end), M);
<<<<<<< HEAD
    fy=polyval(Pi,dx);

    figure
    plot(x,y)
%     hold on
%     plot(dx, fy)
=======
    fy= polyval(Pi,dx);
    
    
>>>>>>> ee520c7bf0def907087a5587872453e16578686d
    
    for i=1:length(Pi)-1
        Pi_d(i)= (N-i)*Pi(i);
    end 
    
    yd= polyval(Pi_d,dx);
    
    
    for j=1:length(Pi_d)-1
        ddpi(j)= (N-j)*Pi_d(j);
    end 
    
    ydd= polyval(ddpi,dx);
    
    g= sqrt(1+(yd).^2);
    
    h= (dx(2)-dx(1));
    
    for i=2:M
        Ltr(i) = (h/2)*( g(1) + g(i) + 2*sum(g(2:1:i-1)));
    end
    
    %calcul de l'erreur 
    dg= (yd.*ydd./g);
    
    E=(h^2/12)*(dg(end)-dg(1));
    
    % calcul de O et vitesse réel 
    
    d = v*ts;
    O = round(Ltr(end)/(d));
    Vr = Ltr(end)/(ts*(O));
    tt = Ltr(end)/Vr;
    dl= Ltr(end)/O;
    % Newton-Raphson
    
    disp(O)
    fn=1000;
    an=x(1);
    bn = an + dl;
    X(1)=an;
    it=0;
    tol=1e-08;
    for i=1:O
        while abs(fn) > tol && it < 1000
            xn = linspace(an, bn, M);
            yd= polyval(Pi_d,xn);
            g= sqrt(1+(yd).^2);
            dn= g(end);
            h1=(xn(2)-xn(1));
            fn = ((g(1) + g(M) + 4*sum(g(2:2:M-1)) + 2*sum(g(3:2:M-1)))*h1/3)-dl;
            bn = bn - ((fn)/dn);
            it=it+1;
            %disp("an : " + an + ", bn : " + bn)
        end 
        it=0;
        X(i)=bn;
        an = bn;
        fn=0.01;
        %disp('hit')
        
    end
    Y = polyval(Pi,X);
    
    disp(X)
    figure
    plot(X, Y)
    hold on
    scatter(dx,fy, "Marker","x")
    legend('N-R', 'Données')
    
end 
 


