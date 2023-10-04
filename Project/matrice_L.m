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
syms Fz(g, Fa, Fb, Fc, ms, mp)
syms Va_eq Vb_eq Vc_eq
syms R Y_d X_d Y_e X_e Y_f X_f


Fa = ((ia_eq^2+b_e1*abs(ia_eq))*sign(ia_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3);
Fb = ((ib_eq^2+b_e1*abs(ib_eq))*sign(ib_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3);
Fc = ((ic_eq^2+b_e1*abs(ic_eq))*sign(ic_eq))/(a_e0+a_e1*z_eq+a_e2*z_eq^2+a_e3*z_eq^3)-1/(a_s0+a_s1*z_eq+a_s2*z_eq^2+a_s3*z_eq^3);

Mx = Fc*XC-Fb*XB+(ms*g*sqrt(Px_eq^2+Py_eq^2)*cos(atan2(Py_eq,Px_eq)));
My = Fa*XA-Fc*XC-Fb*XB+(ms*g*sqrt(Px_eq^2+Py_eq^2)*sin(atan2(Py_eq,Px_eq)));
Fz = g+((Fa+Fb+Fc)/(ms+mp));
dx =-5*g*Ay_eq/7; 
dy =-5*g*Ax_eq/7;
dia= Va_eq-R*ia_eq;
dib= Vb_eq-R*ib_eq;
dic= Vc_eq-R*ic_eq;

PP = [diff(Mx,Ax_eq) diff(Mx,Ay_eq) diff(Mx,z_eq);
      diff(My,Ax_eq) diff(My,Ay_eq) diff(My,z_eq);
      diff(Fz,Ax_eq) diff(Fz,Ay_eq) diff(Fz,z_eq)];
  
SP = [diff(dx,Ax_eq) diff(dx,Ay_eq) diff(dx,z_eq);
      diff(dy,Ax_eq) diff(dy,Ay_eq) diff(dy,z_eq)];
  
CC = [diff(dia,ia_eq) diff(dia,ib_eq) diff(dia,ic_eq);
      diff(dib,ia_eq) diff(dib,ib_eq) diff(dib,ic_eq);
      diff(dic,ia_eq) diff(dic,ib_eq) diff(dic,ic_eq)];


PS = [(diff(Mx,Px_eq)) (diff(Mx,Py_eq)) ; 
      (diff(My,Px_eq)) (diff(My,Py_eq)) ;
      (diff(Fz,Px_eq)) (diff(Fz,Py_eq))];

CV = [(diff(dia,Va_eq)) (diff(dia,Vb_eq)) (diff(dia,Vc_eq)) ; 
      (diff(dib,Va_eq)) (diff(dib,Vb_eq)) (diff(dib,Vc_eq)) ;
      (diff(dic,Va_eq)) (diff(dic,Vb_eq)) (diff(dic,Vc_eq))];

PC = [(diff(Mx,ia_eq)) (diff(Mx,ib_eq)) (diff(Mx,ic_eq)) ; 
      (diff(My,ia_eq)) (diff(My,ib_eq)) (diff(My,ic_eq)) ;
      (diff(Fz,ia_eq)) (diff(Fz,ib_eq)) (diff(Fz,ic_eq))];
  
Tdef=[Y_d -X_d 1;
      Y_e -X_e 1;
      Y_f -X_f 1];
  