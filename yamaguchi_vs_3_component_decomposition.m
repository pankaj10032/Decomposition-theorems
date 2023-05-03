 clc;
clear all;
close all;




% T12=complex(T12_real,T12_imag);
% T13=complex(T13_real,T13_imag);
% T23=complex(T23_real,T23_imag);
% T21=complex(T12_real, -T12_imag);
% T31=complex(T13_real, -T13_imag);
% T32=complex(T23_real,-T23_imag);
%%%%% 20th example%%%%%%
% T11=1.89;
% T12=1.68+1i*0.27;
% T13=-0.35-1i*0.13;
% T21=1.68-1i*0.27;
% T22=4.51;
% T23=-0.46-1i*0.21;
% T31=-0.35+1i*0.13;
% T32=-0.46+1i*0.21;
% T33=0.28;

%%%%%%%% 21th example%%%%%%%
% T11=4.56;
% T12=2.28+1i*0.72;
% T13=0.02+1i*0.67;
% T21=2.28-1i*0.72;
% T22=6.06;
% T23=1.90+1i*0.27;
% T31=0.02*1i*0.67;
% T32=1.90-1i*0.27;
% T33=3.50;

%%%%%%%%%22th example%%%%%%%
T11=0.024;
T12=-0.003-1i*0.004;
T13=0.003-1i*0.000;
T21=conj(T12);
T22=0.006;
T23=-0.001+1i*0.001;
T31=conj(T13);
T32=conj(T23);
T33=0.008;


T23_real=-0.001;
T23_imag=0.001;

T_v=1/4*[2 0 0; 0 1 0; 0 0 1];
T_c1=1/2*[0 0 0;0 1 +1i; 0 -1i 1];
T_c2=1.2*[0 0 0; 0 1 -1i; 0 +1i 1];
T=[T11 T12 T13; T21 T22  T23;  T31 T32  T33];
thetafp_old=atand((T11+T22+T33)*(T11-T22-T33)/(T11*(T22+T33)+ (T11+T22+T33)*(T11+T22+T33)));
T(isnan(T))=0;
eig_val = eig(T,T_c1); %it will have 3 eigenvalues Lambda_1,Lambda_2, Lambda_3.
pc1 =2*T23_imag;   %Here the volume scattering power is equal to the minimum of the eigenvalue


% T11r=T11-0.5*pv1;
% T22r=T22-.25*pv1;
% T33r=T33-0.25*pv1;
T11r=T11;
T22r=T22-pc1/2;
T33r=T33-pc1/2;
 
Tr=[T11r T12 T13; T21 T22r T23; T31 T32 T33r];
tic;
%%%%%%%%%%%%%%%%%%%%%%%%% OAC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetafp_new=atand(trace(Tr)*(T11r-T22r-T33r)/(T11r*(T22r+T33r)+ trace(Tr)*trace(Tr)));

