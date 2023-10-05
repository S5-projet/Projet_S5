clc; 
clear all; 
close all; 


syms Wx_eq Wy_eq % phi' theta'
syms Vx_eq Vy_eq Vz_eq % x' y' z'
syms Fa(ia,  a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1) 
syms Fb(ib,  a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1)
syms Fc(ic,  a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1)
syms Mx(Fb, Fc, Px, Py, XC, XB, ms, g, Ax, Ay)
syms My(Fa, Fb, Fc, Px, Py, XC, XB, XA, ms,g)
syms Fz(g, Fa, Fb, Fc, ms, mp)
syms Va Vb Vc L Pz
syms R 

Ax_des= 0;
Ay_des=0;
Px_des=0;
Py_des=0; 
Pz_des=0;
L=115;
R= 3.6;
g = 9.81;
b_e1 = 13.029359254409743;
ms = 0.008;
mp = 0.442;
Yd = 0.080*cosd(30);
Xd = 0.080*sind(30);
Ye = 0;
Xe = 0.080;
Yf = 0.080*cosd(30);
Xf = 0.080*sind(30);
XA = 0.0952;
XB = XA*sind(30);
XC = XA*sind(30);
a_s0 =1 ;
a_s1 =1 ;
a_s2 =1 ;
a_s3 =1 ;
a_e0 =1 ; 
a_e1 =1 ; 
a_e2 =1 ; 
a_e3 =1 ;
YA = 0;
YB = XA*cosd(30);
YC = XA*cosd(30);
za= YA*Ax-XA*Ay+Pz;
zb= YB*Ax-XB*Ay+Pz;
zc= YC*Ax-XC*Ay+Pz;
ze= Ye*Ax-Xe*Ay+Pz;
zd= Yd*Ax-Xd*Ay+Pz;
zf= Yf*Ax-Xf*Ay+Pz;
ia_eq=-0.03722;
ib_eq=-0.03722;
ic_eq=-0.03722;
Va_eq =-0.1339;
Vb_eq =-0.1339;
Vc_eq =-0.1339;



Fa = ((ia^2+b_e1*abs(ia))*sign(ia))/(a_e0+a_e1*za+a_e2*za^2+a_e3*za^3)-1/(a_s0+a_s1*za+a_s2*za^2+a_s3*za^3);
Fb = ((ib^2+b_e1*abs(ib))*sign(ib))/(a_e0+a_e1*zb+a_e2*zb^2+a_e3*zb^3)-1/(a_s0+a_s1*zb+a_s2*zb^2+a_s3*zb^3);
Fc = ((ic^2+b_e1*abs(ic))*sign(ic))/(a_e0+a_e1*zc+a_e2*zc^2+a_e3*zc^3)-1/(a_s0+a_s1*zc+a_s2*zc^2+a_s3*zc^3);

Mx = Fc*XC-Fb*XB+(ms*g*Py);
My = Fa*XA-Fc*XC-Fb*XB+(ms*g*Px);
Fz = g+((Fa+Fb+Fc)/(ms+mp));
dx =(-5*g*Ay)/7; 
dy =(-5*g*Ax)/7;
dia= (Va-R*ia)/L;
dib= (Vb-R*ib)/L;
dic= (Vc-R*ic)/L;


PP = [diff(Mx,Ax) diff(Mx,Ay) diff(Mx,Pz);
      diff(My,Ax) diff(My,Ay) diff(My,Pz);
      diff(Fz,Ax) diff(Fz,Ay) diff(Fz,Pz)];
  
PP_eq = subs(PP,[ia ib ic Va Vb Vc Ax Ay Px Py Pz],[ia_eq ib_eq ic_eq Va_eq Vb_eq Vc_eq Ax_des Ay_des Px_des Py_des Pz_des]);


SP = [diff(dx,Ax) diff(dx,Ay) diff(dx,Pz);
      diff(dy,Ax) diff(dy,Ay) diff(dy,Pz)];
SP_eq = subs(PP,[ia ib ic Va Vb Vc Ax Ay Px Py Pz],[ia_eq ib_eq ic_eq Va_eq Vb_eq Vc_eq Ax_des Ay_des Px_des Py_des Pz_des]); 
  

CC = [diff(dia,ia) diff(dia,ib) diff(dia,ic);
      diff(dib,ia) diff(dib,ib) diff(dib,ic);
      diff(dic,ia) diff(dic,ib) diff(dic,ic)];


PS = [(diff(Mx,Px)) (diff(Mx,Py)) ; 
      (diff(My,Px)) (diff(My,Py)) ;
      (diff(Fz,Px)) (diff(Fz,Py))];

CV = [(diff(dia,Va)) (diff(dia,Vb)) (diff(dia,Vc)) ; 
      (diff(dib,Va)) (diff(dib,Vb)) (diff(dib,Vc)) ;
      (diff(dic,Va)) (diff(dic,Vb)) (diff(dic,Vc))];

PC = [(diff(Mx,ia)) (diff(Mx,ib)) (diff(Mx,ic)) ; 
      (diff(My,ia)) (diff(My,ib)) (diff(My,ic)) ;
      (diff(Fz,ia)) (diff(Fz,ib)) (diff(Fz,ic))];
  
Tdef=[Yd -Xd 1;
      Ye -Xe 1;
      Yf -Xf 1];
  
I=eye(4);
C1=zeros(3,4);
C2= zeros(3,1);
C3=zeros(3,3);
B1=zeros(10,3);

A = sym(zeros(13, 13));
revue3x3 = eye(3);
revue2x2 = eye(2);

% Insérez les petites matrices à des emplacements spécifiques dans la grande matrice
A(4:6, 1:3) = PP; 
A(1:3, 4:6) = revue3x3;
A(4:6, 7:8) = PS; 
A(4:6, 11:13) = PC;
A(7:8, 9:10) = revue2x2; 
A(9:10, 1:3) = SP;
A(11:13, 11:13) = CC;

% Affichez la matrice X dans la console en utilisant la fonction disp
disp(A);


B=[B1
   CV];

C=[Tdef C3 C3 C3 C2;
   0 0 0 0 0 0 1 0 0 0 0 0 0
   0 0 0 0 0 0 0 1 0 0 0 0 0
   0 0 0 0 0 0 0 0 1 0 0 0 0
   0 0 0 0 0 0 0 0 0 1 0 0 0];


  