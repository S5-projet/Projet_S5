clc 
close all 

syms Wx_eq Wy_eq % phi' theta'
syms Vx_eq Vy_eq Vz_eq % x' y' z'
syms Ax Pz Px Py ia_eq ib_eq ic_eq Ay 

% Constante
g = 9.81;  
XA = 0.0952;
XB =- XA*sind(30);
XC =- XA*sind(30);
YA = 0;
YB = XA*cosd(30);
YC = -XA*cosd(30);
Yd = 0.080*cosd(30);
Xd = 0.080*sind(30);
Ye = 0;  
Xe = -0.080;
Yf = -0.080*cosd(30);
Xf = 0.080*sind(30);
b_e1 = 13.029359254409743; 
ms = 0.0;
mp = 0.442;
R= 3.6;

%variable
Jx= 1347/1000^2;
Jy=Jx;
Px= 0;
Py= 0;
Pz=0.011;
Ax=deg2rad(0);
Ay=deg2rad(0.0);
z_eq= 0.011; 
a_s0 =0.0609;
a_s1 =18.1233;
a_s2 =3.1142e+03;
a_s3 =-2.3281e+05;
a_s4 =3.9161e+07;
a_s5 =-1.0305e+09;
a_s6 =-3.1858e+07;
a_e0 =0.000013398315537 * 1.0e+5 ; 
a_e1 =0.003556457401249 * 1.0e+5 ; 
a_e2 =0.018168930792891 * 1.0e+5 ; 
a_e3 =7.043847497815494 * 1.0e+5 ;
za = YA*Ax-XA*Ay+Pz
zb = YB*Ax-XB*Ay+Pz
zc = YC*Ax-XC*Ay+Pz
% ze= Ye*Ax-Xe*Ay+Pz;
% zd= Yd*Ax-Xd*Ay+Pz;
% zf= Yf*Ax-Xf*Ay+Pz;

% equation pour trouver ia, ib, ic
Fa = ((ia_eq^2+b_e1*abs(ia_eq))*sign(ia_eq))/(a_e0+a_e1*za+a_e2*za^2+a_e3*za^3)-(1/(a_s0+a_s1*za+a_s2*za^2+a_s3*za^3+a_s4*za^4+a_s5*za^5+a_s6*za^6));
Fb = ((ib_eq^2+b_e1*abs(ib_eq))*sign(ib_eq))/(a_e0+a_e1*zb+a_e2*zb^2+a_e3*zb^3)-(1/(a_s0+a_s1*zb+a_s2*zb^2+a_s3*zb^3+a_s4*zb^4+a_s5*zb^5+a_s6*zb^6));
Fc = ((ic_eq^2+b_e1*abs(ic_eq))*sign(ic_eq))/(a_e0+a_e1*zc+a_e2*zc^2+a_e3*zc^3)-(1/(a_s0+a_s1*zc+a_s2*zc^2+a_s3*zc^3+a_s4*zc^4+a_s5*zc^5+a_s6*zc^63));




My = ((-Fa*XA-Fc*XC-Fb*XB+(ms*g*Px)))/Jx ==0;
Mx = ((Fc*YC+Fb*YB-(ms*g*Py)))/Jy ==0;
Fz =(Fa + Fb + Fc + mp*g + ms*g) / mp ==0;

[ia_eq, ib_eq, ic_eq] = vpasolve(Mx, My, Fz, ia_eq, ib_eq, ic_eq)

Va_eq = R*ia_eq;
Vb_eq = R*ib_eq;
Vc_eq = R*ic_eq;

Va_eq_simu = double(Va_eq);%Besoin pour simulink (on peux pas faire de simulation avec des variables symboliques)
Vb_eq_simu = double(Vb_eq);%"
Vc_eq_simu = double(Vc_eq);%"

Fea=((ia_eq^2+b_e1*abs(ia_eq))*sign(ia_eq))/(a_e0+a_e1*za+a_e2*za^2+a_e3*za^3)

Fsa=-(1/(a_s0+a_s1*za+a_s2*za^2+a_s3*za^3+a_s4*za^4+a_s5*za^5+a_s6*za^6))

Feb = ((ib_eq^2+b_e1*abs(ib_eq))*sign(ib_eq))/(a_e0+a_e1*zb+a_e2*zb^2+a_e3*zb^3)
Fsb = (1/(a_s0+a_s1*zb+a_s2*zb^2+a_s3*zb^3+a_s4*zb^4+a_s5*zb^5+a_s6*zb^6))

Fec = ((ic_eq^2+b_e1*abs(ic_eq))*sign(ic_eq))/(a_e0+a_e1*zc+a_e2*zc^2+a_e3*zc^3)
Fsc = (1/(a_s0+a_s1*zc+a_s2*zc^2+a_s3*zc^3+a_s4*zc^4+a_s5*zc^5+a_s6*zc^63))

x_ini= Px
y_ini=Py
z_ini=Pz
v_ini_phi=0
v_ini_theta=0
v_ini_z=0
xs_ini=0
ys_ini=0
vxs_ini=0
vys_ini=0

% x_ini= Px
% y_ini=Py
% z_ini=Pz
% v_ini_phi=0
% v_ini_theta=0
% v_ini_z=0
% xs_ini=0.003
% ys_ini=0.050
% vxs_ini=0
% vys_ini=0








