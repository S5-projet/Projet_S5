close all
clc



syms Wx_eq Wy_eq % phi' theta'
syms Vx_eq Vy_eq Vz_eq % x' y' z'
syms Fa(ia,  a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1) 
syms Fb(ib,  a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1)
syms Fc(ic,  a_s0, a_s1, a_s2, a_s3, a_e0, a_e1, a_e2, a_e3, b_e1)
syms Mx(Fb, Fc, Px, Py, XC, XB, ms, g, Ax, Ay)
syms My(Fa, Fb, Fc, Px, Py, XC, XB, XA, ms,  g)
syms Fz(g, Fa, Fb, Fc, ms, mp)
syms Va Vb Vc L Pz I_phi I_theta I_z V_phi V_theta V_z
syms R L

Ax_des=0;
Ay_des=0;
Px_des=0;
Py_des=0; 
Pz_des=0.015;
L=0.115;
R= 3.6;
g = 9.81;
b_e1 = 13.029359254409743;
ms = 0.008;
mp = 0.442;          
rs = 0.0039;           
Js = 2*ms*rs^2/5; 
Jx= 1347/1000^2;
Jy=Jx;
Yd = 0.080*cosd(30);
Xd = 0.080*sind(30);
Ye = 0;  
Xe = -0.080;
Yf = -0.080*cosd(30);
Xf = 0.080*sind(30);
XA = 0.0952;
XB = -XA*sind(30);
XC = -XA*sind(30);
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
YA = 0;
YB = XA*cosd(30);
YC =-XA*cosd(30);
za= YA*Ax-XA*Ay+Pz;
zb= YB*Ax-XB*Ay+Pz;
zc= YC*Ax-XC*Ay+Pz;
ze= Ye*Ax-Xe*Ay+Pz;
zd= Yd*Ax-Xd*Ay+Pz;
zf= Yf*Ax-Xf*Ay+Pz;
ia_eq= -0.5440584309597635796798512805225;
ib_eq= -0.5440584309597635796798512805225;
ic_eq= -0.5440584309597635796798512805225;
Va_eq= -1.958610351455148886847464609881;
Vb_eq= -1.958610351455148886847464609881;
Vc_eq= -1.958610351455148886847464609881;



Fa = ((ia^2+b_e1*abs(ia))*sign(ia))/(a_e0+a_e1*za+a_e2*za^2+a_e3*za^3)-(1/(a_s0+a_s1*za+a_s2*za^2+a_s3*za^3+a_s4*za^4+a_s5*za^5+a_s6*za^6));
Fb = ((ib^2+b_e1*abs(ib))*sign(ib))/(a_e0+a_e1*zb+a_e2*zb^2+a_e3*zb^3)-(1/(a_s0+a_s1*zb+a_s2*zb^2+a_s3*zb^3+a_s4*zb^4+a_s5*zb^5+a_s6*zb^6));
Fc = ((ic^2+b_e1*abs(ic))*sign(ic))/(a_e0+a_e1*zc+a_e2*zc^2+a_e3*zc^3)-(1/(a_s0+a_s1*zc+a_s2*zc^2+a_s3*zc^3+a_s4*zc^4+a_s5*zc^5+a_s6*zc^6));

My = ((-Fa*XA-Fc*XC-Fb*XB+(ms*g*Px))*cosd(Ay))/Jx;
Mx = ((Fc*YC+Fb*YB-(ms*g*Py))*cosd(Ax))/Jy;
Fz = g +((Fa+Fb+Fc)/(ms+mp));
dx = -(ms*g*Ay) / (ms + Js/rs^2);
dy = (ms*g*Ax) / (ms + Js/rs^2);
dia= (Va-R*ia)/L;
dib= (Vb-R*ib)/L;
dic= (Vc-R*ic)/L;


PP = [diff(Mx,Ax) diff(Mx,Ay) diff(Mx,Pz);
      diff(My,Ax) diff(My,Ay) diff(My,Pz);
      diff(Fz,Ax) diff(Fz,Ay) diff(Fz,Pz)];
  
PP_eq = double(subs(PP,[ia ib ic Va Vb Vc Ax Ay Px Py Pz],[ia_eq ib_eq ic_eq Va_eq Vb_eq Vc_eq Ax_des Ay_des Px_des Py_des Pz_des]));


SP = [diff(dx,Ax) diff(dx,Ay) diff(dx,Pz);
      diff(dy,Ax) diff(dy,Ay) diff(dy,Pz)];

SP_eq = double(subs(SP,[ia ib ic Va Vb Vc Ax Ay Px Py Pz],[ia_eq ib_eq ic_eq Va_eq Vb_eq Vc_eq Ax_des Ay_des Px_des Py_des Pz_des])); 
  

CC = [diff(dia,ia) diff(dia,ib) diff(dia,ic) ;
      diff(dib,ia) diff(dib,ib) diff(dib,ic) ;
      diff(dic,ia) diff(dic,ib) diff(dic,ic)];

CC_eq = double(subs(CC,[ia ib ic Va Vb Vc Ax Ay Px Py Pz],[ia_eq ib_eq ic_eq Va_eq Vb_eq Vc_eq Ax_des Ay_des Px_des Py_des Pz_des])); 

PS = [(diff(Mx,Px)) (diff(Mx,Py)) ; 
      (diff(My,Px)) (diff(My,Py)) ;
      (diff(Fz,Px)) (diff(Fz,Py))];



CV = [(diff(dia,Va)) (diff(dia,Vb)) (diff(dia,Vc)) ; 
      (diff(dib,Va)) (diff(dib,Vb)) (diff(dib,Vc)) ;
      (diff(dic,Va)) (diff(dic,Vb)) (diff(dic,Vc))];

CV_eq = double(subs(CV,[ia ib ic Va Vb Vc Ax Ay Px Py Pz],[ia_eq ib_eq ic_eq Va_eq Vb_eq Vc_eq Ax_des Ay_des Px_des Py_des Pz_des])); 

PC = [(diff(Mx,ia)) (diff(Mx,ib)) (diff(Mx,ic)) ; 
      (diff(My,ia)) (diff(My,ib)) (diff(My,ic)) ;
      (diff(Fz,ia)) (diff(Fz,ib)) (diff(Fz,ic))];

PC_eq = double(subs(PC,[ia ib ic Va Vb Vc Ax Ay Px Py Pz],[ia_eq ib_eq ic_eq Va_eq Vb_eq Vc_eq Ax_des Ay_des Px_des Py_des Pz_des]));  



Tdef=[Yd -Xd 1;
      Ye -Xe 1;
      Yf -Xf 1];
  
I=eye(4);
C1=zeros(3,4);
C2=zeros(3,1);
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
%disp(A);


B=[B1
   CV];

C=[Tdef C3 C3 C3 C2;
   0 0 0 0 0 0 1 0 0 0 0 0 0
   0 0 0 0 0 0 0 1 0 0 0 0 0
   0 0 0 0 0 0 0 0 1 0 0 0 0
   0 0 0 0 0 0 0 0 0 1 0 0 0];

D = [0];


U = [YA YB YC;
    -XA -XB -XC;
     1 1 1];
 
Uinv= inv(U);

I = [I_phi; I_theta; I_z];
V = [V_phi; V_theta; V_z];
V_new= Uinv*V;
I_new= Uinv*I; 
ia_new= I_new(1);
ib_new= I_new(2); 
ic_new= I_new(3);
Va_new= V_new(1);
Vb_new= V_new(2); 
Vc_new= V_new(3);

Fa_new = subs(Fa,[ia],[ia_new]);
Fb_new = subs(Fb,[ib],[ib_new]);
Fc_new = subs(Fc,[ic],[ic_new]);

My_new = ((Fa_new*XA-Fc_new*XC-Fb_new*XB)*cosd(Ay));
Mx_new = ((Fc_new*XC-Fb_new*XB)*cosd(Ax));
Fz_new = g+((Fa_new+Fb_new+Fc_new)/(mp));
dia_new= (Va_new-R*ia_new)/L;
dib_new= (Vb_new-R*ib_new)/L;
dic_new= (Vc_new-R*ic_new)/L;


I_ptz = U*[ia_eq ib_eq ic_eq]';  
V_ptz = U*[Va_eq Vb_eq Vc_eq]';
IP_eq = I_ptz(1);
It_eq = I_ptz(2);
Iz_eq = I_ptz(3);
Vp_eq = V_ptz(1);
Vt_eq = V_ptz(2);
Vz_eq = V_ptz(3);



PP_new = [diff(Mx_new,Ax) diff(Mx_new,Ay) diff(Mx_new,Pz);
          diff(My_new,Ax) diff(My_new,Ay) diff(My_new,Pz);
          diff(Fz_new,Ax) diff(Fz_new,Ay) diff(Fz_new,Pz)];

PP_eq_new = double(subs(PP_new,[I_phi I_theta I_z V_phi V_theta V_z Ax Ay Px Py Pz],[IP_eq It_eq Iz_eq Vp_eq Vt_eq Vz_eq Ax_des Ay_des Px_des Py_des Pz_des]));
      
      
PC_new = [(diff(Mx_new,I_phi)) (diff(Mx_new,I_theta)) (diff(Mx_new,I_z)) ; 
          (diff(My_new,I_phi)) (diff(My_new,I_theta)) (diff(My_new,I_z)) ;
          (diff(Fz_new,I_phi)) (diff(Fz_new,I_theta)) (diff(Fz_new,I_z))];
      
PC_eq_new = double(subs(PC_new,[I_phi I_theta I_z V_phi V_theta V_z Ax Ay Px Py Pz],[IP_eq It_eq Iz_eq Vp_eq Vt_eq Vz_eq Ax_des Ay_des Px_des Py_des Pz_des]));
            
          
CC_new = [diff(dia_new,I_phi) diff(dia_new,I_theta) diff(dia_new,I_z) ;
          diff(dib_new,I_phi) diff(dib_new,I_theta) diff(dib_new,I_z) ;
          diff(dic_new,I_phi) diff(dic_new,I_theta) diff(dic_new,I_z)];
      
CC_eq_new = double(subs(CC_new,[I_phi I_theta I_z V_phi V_theta V_z Ax Ay Px Py Pz],[IP_eq It_eq Iz_eq Vp_eq Vt_eq Vz_eq Ax_des Ay_des Px_des Py_des Pz_des]));
           
      
CV_new = [(diff(dia_new,V_phi)) (diff(dia_new,V_theta)) (diff(dia_new,V_z)) ; 
          (diff(dib_new,V_phi)) (diff(dib_new,V_theta)) (diff(dib_new,V_z)) ;
          (diff(dic_new,V_phi)) (diff(dic_new,V_theta)) (diff(dic_new,V_z))];
      
CV_eq_new = double(subs(CV_new,[I_phi I_theta I_z V_phi V_theta V_z Ax Ay Px Py Pz],[IP_eq It_eq Iz_eq Vp_eq Vt_eq Vz_eq Ax_des Ay_des Px_des Py_des Pz_des]));
      


% valeur équilibre 

zd_eq=0.015;
ze_eq=0.015;
zf_eq=0.015;
Px_eq=0;
Py_eq=0; 
Vx_eq=0;
Vy_eq=0;

y_eq=[zd_eq ze_eq zf_eq Px_eq Py_eq Vx_eq Vy_eq]';
