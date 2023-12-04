clc 
close all 
clear all
x=[0 1 2 3]';
y=[0 1 2 3]';
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

M=101;
dx = linspace(x(1), x(end), M);
fy=polyval(Pi,dx);



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
Vr = Ltr(end)*(ts/(O));
dl= Ltr(end)/O;
% Newton-Raphson


fn(1)=0;
an=x(1);
X(1)=an;
for i=1:O
    bn = an + 0.001;
    xn = linspace(an, bn, M);
    yd= polyval(Pi_d,xn);
    g= sqrt(1+(yd).^2);
    dn= g;
    h1=(xn(2)-xn(1));
    fn(((i-1)/2)+1)=(g(1)+g(i)+4*sum(g(2:2:i-1))+2*sum(g(3:2:i-1)))*h/3;
    it=0;
    tol=1e-08;


while abs(fn)> tol
    bn = bn-((fn-dl)./dn);
    xn = linspace(an, bn, M);
    yd= polyval(Pi_d,xn);
    g= sqrt(1+(yd).^2);
    dn= g;
    h1=(xn(2)-xn(1));
    fn(((i-1)/2)+1)=(g(1)+g(i)+4*sum(g(2:2:i-1))+2*sum(g(3:2:i-1)))*h/3;
    it=it+1;
end 
X(i)=bn;
an = bn;
end
Y = polyval(flip(Pi),X)
end 
 
% dl = racine xn^2+ y(xn)^2

