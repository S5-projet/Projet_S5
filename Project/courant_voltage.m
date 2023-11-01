syms Wx_eq Wy_eq % phi' theta'
syms Vx_eq Vy_eq Vz_eq % x' y' z'
syms Fa(ia_eq, z_eq, a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1) 
syms Fb(ib_eq, z_eq, a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1)
syms Fc(ic_eq, z_eq, a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1)
syms Mx(Fb, Fc, Px_eq, Py_eq, XC, XB, ms, g, Ax_eq, Ay_eq)
syms My(Fa, Fb, Fc, Px_eq, Py_eq, XC, XB, XA, ms,g)
syms fz(g, Fa, Fb, Fc, ms, mp)

% Constante
g = 9.81;  
XA = 0.0952;
XB = XA*sind(30);
XC = XA*sind(30);
YA = 0;
YB = XA*cosd(30);
YC = XA*cosd(30);
b_e1 = 13.029359254409743; 
ms = 0.008;
mp = 0.442;
R= 3.6;

%variable
Px_eq= 0;
Py_eq= 0;
z_eq= 0.010; 
a_s0 =0.000000000001284 * 1.0e+11 ;
a_s1 =-0.000000000792998 * 1.0e+11;
a_s2 =0.000000381718748 * 1.0e+11 ;
a_s3 =-0.000052597553068 * 1.0e+11;
a_s4 =0.003775582389408 * 1.0e+11 ;
a_s5 =-0.122041077197557 * 1.0e+11;
a_s6 =1.535057200314033 * 1.0e+11 ;
a_e0 =0.000013398315537 * 1.0e+5 ; 
a_e1 =0.003556457401249 * 1.0e+5 ; 
a_e2 =0.018168930792891 * 1.0e+5 ; 
a_e3 =7.043847497815494 * 1.0e+5 ;

% equation pour trouver ia, ib, ic
+Fa = ((ia_eq^2+b_e1*abs(ia_eq))*sign(ia_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3+a_s4*z_eq^4+a_s5*z_eq^5+a_s6*z_eq^6);
Fb = ((ib_eq^2+b_e1*abs(ib_eq))*sign(ib_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3+a_s4*z_eq^4+a_s5*z_eq^5+a_s6*z_eq^6);
Fc = ((ic_eq^2+b_e1*abs(ic_eq))*sign(ic_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3+a_s4*z_eq^4+a_s5*z_eq^5+a_s6*z_eq^63);

Mx = Fc*XC-Fb*XB+(ms*g*Py_eq)==0;
My = Fa*XA-Fc*XC-Fb*XB+(ms*g*Px_eq)==0;
Fz= g+((Fa+Fb+Fc)/(ms+mp))==0;

[ia_eq, ib_eq, ic_eq] = vpasolve(Mx, My, Fz, ia_eq, ib_eq, ic_eq)

Va_eq = R*ia_eq;
Vb_eq = R*ib_eq;
Vc_eq = R*ic_eq;










