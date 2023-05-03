 clc;
clear all;
close all;



% 
% T12=complex(T12_real,T12_imag);
% T13=complex(T13_real,T13_imag);
% T23=complex(T23_real,T23_imag);
% T21=complex(T12_real, -T12_imag);
% T31=complex(T13_real, -T13_imag);
% T32=complex(T23_real,-T23_imag);
%%%% 20th example%%%%%%
% T11=1.89;
% T12=1.68+1i*0.27;
% T13=-0.35-1i*0.13;
% T21=1.68-1i*0.27;
% T22=4.51;
% T23=-0.46-1i*0.21;
% T31=-0.35+1i*0.13;
% T32=-0.46+1i*0.21;
% T33=0.28;

%%%%%%% 21th example%%%%%%%
% T11=4.56;
% T12=2.28+1i*0.72;
% T13=0.02+1i*0.67;
% T21=2.28-1i*0.72;
% T22=6.06;
% T23=1.90+1i*0.27;
% T31=0.02*1i*0.67;
% T32=1.90-1i*0.27;
% T33=3.50;

% %%%%%%%%%22th example%%%%%%%
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

T=[T11 T12 T13; T21 T22 T23; T31 T32 T33];
span_T=trace(T);
trace_T=trace(T);
theta_FP_T = real(atand((span_T*(T(1,1) - T(2,2) - T(3,3)))/((T(1,1)*(T(2,2) + T(3,3)))+((span_T^2)))));

theta=1/4*atand((2*T23_imag)/(T22-T33));
R_x=[1 0 0; 0 cos(theta) 1i*sin(theta); 0 1i*sin(theta) cos(theta)];
T_x=R_x*T*R_x';
span_T_x=trace(T_x);
trace_T_x=trace(T_x);
theta_FP_T_x = real(atand((span_T_x*(T_x(1,1) - T_x(2,2) - T_x(3,3)))/((T_x(1,1)*(T_x(2,2) + T_x(3,3)))+((span_T_x^2)))));
ps=trace_T*(1+sind(2*theta_FP_T))/2;
pd=trace_T*(1-sind(2*theta_FP_T))/2;
