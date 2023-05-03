clc;
clear all;
close all;




T11=4.56;
T12=2.28+1i*0.72;
T13=0.02+1i*0.67;
T21=2.28-1i*0.72;
T22=6.06;
T23=1.90+1i*0.27;
T31=0.02-1i*0.72;
T32=1.90-1i*0.72;
T33=3.50;

T_v=1/4*[2 0 0; 0 1 0; 0 0 1];
T23_real=1.90;
T23_imag=0.27;
%T=[T11 complex(T12_real,T12_imag) complex(T13_real,T13_imag);...
  %conj(complex(T12_real,T12_imag)) T22 complex(T23_real,T23_imag);...
  %conj(complex(T13_real,T13_imag)) conj(complex(T23_real,T23_imag)) T33];
T=[T11 T12 T13; T21 T22  T23;  T31 T32  T33];
T(isnan(T))=0;
eig_val = eig(T,T_v); %it will have 3 eigenvalues Lambda_1,Lambda_2, Lambda_3.
pv1 = min(eig_val);  %Here the volume scattering power is equal to the minimum of the eigenvalue


T11r=T11;
T22r=T22;
T33r=T33;

Tr=[T11r T12 T13; T21 T22r T23; T31 T32 T33r];
tic;
%%%%%%%%%%%%%%%%%%%%%%%%% OAC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1=(0.25)*atan(2*T23_real./(T22r-T33r));
thetafp1=atand((T11r+T22r+T33r)*(T11r-T22r-T33r)/(T11r*(T22r+T33r)+ (T11r+T22r+T33r)*(T11r+T22r+T33r)));
ps1=(T11r+T22r+T33r)*(1+sind(2*thetafp1))/2;
pd1=(T11r+T22r+T33r)*(1-sind(2*thetafp1))/2;
T11n=T11r;
T12n=T12.*cos(2*v1)+T13.*sin(2*v1);
T13n=-T12.*sin(2*v1)+T13.*cos(2*v1);
T23n=1/2*(T33r-T22r).*sin(4*v1)+(T23_real).*cos(4*v1)+1i*(T23_imag);
T22n=T22r.*(cos(2*v1).*cos(2*v1)) + T33r.*(sin(2*v1).*sin(2*v1))+(T23_real.*sin(4*v1));
T33n=T22r.*(sin(2*v1).*sin(2*v1)) + T33r.*(cos(2*v1).*cos(2*v1))- (T23_real.*sin(4*v1));
thetafp_middle=atand((T11n+T22n+T33n)*(T11n-T22n-T33n)/(T11n*(T22n+T33n)+ (T11n+T22n+T33n)*(T11n+T22n+T33n)));
ps_middle=(T11n+T22n+T33n)*(1+sind(2*thetafp_middle))/2;
pd_middle=(T11n+T22n+T33n)*(1-sind(2*thetafp_middle))/2;