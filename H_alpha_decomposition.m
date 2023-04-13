clc;

clear all;
close all;

fid1 = fopen("C:\Alos_San_Ship1\T3\T11.bin");
fid2 = fopen("C:\Alos_San_Ship1\T3\T22.bin");
fid3 = fopen("C:\Alos_San_Ship1\T3\T33.bin");
fid4 = fopen("C:\Alos_San_Ship1\T3\T12_real.bin");
fid5 = fopen("C:\Alos_San_Ship1\T3\T12_imag.bin");
fid6 = fopen("C:\Alos_San_Ship1\T3\T13_real.bin");
fid7 = fopen("C:\Alos_San_Ship1\T3\T13_imag.bin");
fid8 = fopen("C:\Alos_San_Ship1\T3\T23_real.bin");
fid9 = fopen("C:\Alos_San_Ship1\T3\T23_imag.bin");

Nligfin=293;
Ncol=303;


for lig = 1 : Nligfin
    T11(lig,1:Ncol) = fread(fid1, Ncol, 'float32');
    T22(lig,1:Ncol) = fread(fid2, Ncol, 'float32');
    T33(lig,1:Ncol) = fread(fid3, Ncol, 'float32');
    T12_real(lig,1:Ncol) = fread(fid4, Ncol, 'float32');
    T12_imag(lig,1:Ncol) = fread(fid5, Ncol, 'float32');
    T13_real(lig,1:Ncol) = fread(fid6, Ncol, 'float32');
    T13_imag(lig,1:Ncol) = fread(fid7, Ncol, 'float32');
    T23_real(lig,1:Ncol) = fread(fid8, Ncol, 'float32');
    T23_imag(lig,1:Ncol) = fread(fid9, Ncol, 'float32');
end


span=T11+T22+T33;
T12=complex(T12_real,T12_imag);
T13=complex(T13_real,T13_imag);
T23=complex(T23_real,T23_imag);
T21=complex(T12_real, -T12_imag);
T31=complex(T13_real, -T13_imag);
T32=complex(T23_real,-T23_imag);

for i=1:Nligfin
    for j=1:Ncol   
T=[T11(i,j) T12(i,j) T13(i,j); T21(i,j) T22(i,j)  T23(i,j);  T31(i,j) T32(i,j)  T33(i,j)];
T(isnan(T))=0;
[V,D] = eig(T);
lambda_1=D(3,3);
lambda_2=D(2,2);
lambda_3=D(1,1);

%pseudo probabilities
p1=lambda_1/(lambda_1+lambda_2+lambda_3);
p2=lambda_2/(lambda_1+lambda_2+lambda_3);
p3=lambda_3/(lambda_1+lambda_2+lambda_3);

%entropy calculation
H(i,j) =-p1*log(p1)/log(3) -p2*log(p2)/log(3)-p3*log(p3)/log(3);

%Here our every element of eigenvector U have both real and complex parts
%V(1,1)=a+jb
%cos(alpha_1)=sqrt(a*a+b*b)

real_11=real(V(1,1));
imag_11=imag(V(1,1));
alpha_3=acosd(sqrt(real_11*real_11 + imag_11*imag_11));

real_12=real(V(1,2));
imag_12=imag(V(1,2));
alpha_2=acosd(sqrt(real_12*real_12 + imag_12*imag_12));

real_13=real(V(1,3));
imag_13=imag(V(1,3));
alpha_1=acosd(sqrt(real_13*real_13 + imag_13*imag_13));

%alpha calculation
alpha(i,j)=p1*alpha_1+p2*alpha_2+p3*alpha_3;
    end
end


%now we have to classifiy our surface on the basis of alpha and entropy in
%the 9 diffrent parts

for i=1:Nligfin
    for j=1:Ncol
            if alpha(i,j)<42.5 && H(i,j)<0.5 %bragg surface/very smooth land surfaces
                psm1(i,j,1)=0;
                psm1(i,j,2)=0;
                psm1(i,j,3)=1;
            elseif alpha(i,j)>42.5 && alpha(i,j)<47.5 && H(i,j)<0.5
                psm1(i,j,1)=0;
                psm1(i,j,2)=0;
                psm1(i,j,3)=0;
            elseif  alpha(i,j)> 47.5 && H(i,j)<0.5
                psm1(i,j,1)=1;
                psm1(i,j,2)=1;
                psm1(i,j,3)=1;
            elseif alpha(i,j) < 40 && H(i,j)<0.9 && H(i,j)>0.5
                psm1(i,j,1)=1;
                psm1(i,j,2)=0;
                psm1(i,j,3)=0;
            elseif alpha(i,j)< 40 &&  alpha(i,j)> 50 && H(i,j)<0.9 && H(i,j)>0.5
                psm1(i,j,1)=0;
                psm1(i,j,2)=1;
                psm1(i,j,3)=0;
            elseif alpha(i,j)> 50 &&  H(i,j)<0.9 && H(i,j)>0.5
                psm1(i,j,1)=1;
                psm1(i,j,2)=1;
                psm1(i,j,3)=0;
            elseif alpha(i,j) < 40 && H(i,j)>0.9
                psm1(i,j,1)=1;
                psm1(i,j,2)=0;
                psm1(i,j,3)=1;
            elseif alpha(i,j)> 40 && alpha(i,j)< 55 && H(i,j)>0.9
                psm1(i,j,1)=0;
                psm1(i,j,2)=1;
                psm1(i,j,3)=1;
            elseif alpha(i,j)>55 && H(i,j)>0.9
                psm1(i,j,1)=0.25;
                psm1(i,j,2)=0.25;
                psm1(i,j,3)=0;
            end
    end
end



