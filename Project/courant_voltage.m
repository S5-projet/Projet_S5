clc; 
clear all; 
close all; 


syms Wx_eq Wy_eq % phi' theta'
syms Vx_eq Vy_eq Vz_eq % x' y' z'
syms Fa(ia_eq, z_eq, a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1) 
syms Fb(ib_eq, z_eq, a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1)
syms Fc(ic_eq, z_eq, a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1)
syms Mx(Fb, Fc, Px_eq, Py_eq, XC, XB, ms, g, Ax_eq, Ay_eq)
syms My(Fa, Fb, Fc, Px_eq, Py_eq, XC, XB, XA, ms,g)
syms fz(g, Fa, Fb, Fc, ms, mp)

%variable qui ne change pas 
g = 9.81;  
XA = 0.0952;
XB = XA*sind(30);
XC = XA*sind(30);
b_e1 = 13.029359254409743; 
ms = 0.008;
mp = 0.442;
R= 3.6;

%variable qui peux changer
Px_eq= 0;
Py_eq= 0;
z_eq= 0.030; 
a_s0 =1 ;
a_s1 =1 ;
a_s2 =1 ;
a_s3 =1 ;
a_e0 =1 ; 
a_e1 =1 ; 
a_e2 =1 ; 
a_e3 =1 ;

% equation pour trouver ia, ib, ic
Fa = ((ia_eq^2+b_e1*abs(ia_eq))*sign(ia_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3);
Fb = ((ib_eq^2+b_e1*abs(ib_eq))*sign(ib_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3);
Fc = ((ic_eq^2+b_e1*abs(ic_eq))*sign(ic_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3);

Mx = Fc*XC-Fb*XB+(ms*g*sqrt(Px_eq^2+Py_eq^2)*cos(atan2(Py_eq,Px_eq)))==0;
My = Fa*XA-Fc*XC-Fb*XB+(ms*g*sqrt(Px_eq^2+Py_eq^2)*sin(atan2(Py_eq,Px_eq)))==0;
Fz= g+((Fa+Fb+Fc)/(ms+mp))==0;

[ia_eq, ib_eq, ic_eq] = vpasolve(Mx, My, Fz, ia_eq, ib_eq, ic_eq);

Va_eq = R*ia_eq;
Vb_eq = R*ib_eq;
Vc_eq = R*ic_eq;










