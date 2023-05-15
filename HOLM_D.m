clc;
clear all;
close all;


fid1 = fopen("D:\Alos_San_Ship1\T3\T11.bin");
fid2 = fopen("D:\Alos_San_Ship1\T3\T22.bin");
fid3 = fopen("D:\Alos_San_Ship1\T3\T33.bin");
fid4 = fopen("D:\Alos_San_Ship1\T3\T12_real.bin");
fid5 = fopen("D:\Alos_San_Ship1\T3\T12_imag.bin");
fid6 = fopen("D:\Alos_San_Ship1\T3\T13_real.bin");
fid7 = fopen("D:\Alos_San_Ship1\T3\T13_imag.bin");
fid8 = fopen("D:\Alos_San_Ship1\T3\T23_real.bin");
fid9 = fopen("D:\Alos_San_Ship1\T3\T23_imag.bin");

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

c=T11-T22;
span=T11+T22+T33;
T12=complex(T12_real,T12_imag);
T13=complex(T13_real,T13_imag);
T23=complex(T23_real,T23_imag);
T21=complex(T12_real, -T12_imag);
T31=complex(T13_real, -T13_imag);
T32=complex(T23_real,-T23_imag);
%T21=conj(T12);
%T31=conj(T11);
%T32=conj(T23);
tic
for i=1:Nligfin
    for j=1:Ncol   
T=[T11(i,j) T12(i,j) T13(i,j); T21(i,j) T22(i,j)  T23(i,j);  T31(i,j) T32(i,j)  T33(i,j)];
[V,D] = eig(T);
Lambda_1=D(3,3);
Lambda_2=D(2,2);
Lambda_3=D(1,1);
pv1(i,j)=3*Lambda_3;
if c(i,j)>0
    ps1(i,j)=Lambda_1-Lambda_2;
    pd1(i,j)=2*(Lambda_2 - Lambda_3);
else         
    ps1(i,j)=2*(Lambda_2 -Lambda_3);
    pd1(i,j)= Lambda_1-Lambda_2;
end
    end
end

toc
psm2=cat(3,pd1,pv1,ps1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nligfin
    for j=1:Ncol
        if ps1(i,j)>pd1(i,j) && ps1(i,j)>pv1(i,j)       %% zone 8 blue%%
           psm1(i,j,1)=0;
           psm1(i,j,2)=0;
           psm1(i,j,3)=1;
        elseif pd1(i,j)>=ps1(i,j) && pd1(i,j)>=pv1(i,j)       %% zone 8 red%%
           psm1(i,j,1)=1;
           psm1(i,j,2)=0;
           psm1(i,j,3)=0;
        elseif pv1(i,j)>ps1(i,j) && pv1(i,j)>pd1(i,j)   %% zone 6 green%%
           psm1(i,j,1)=0;
           psm1(i,j,2)=1;
           psm1(i,j,3)=0;
           end
     end
end

%%%% for detection ships by applying mask%%%%%
for i=1:Nligfin
    for j=1:Ncol
        if ps1(i,j)>pd1(i,j) && ps1(i,j)>pv1(i,j)
             psm1(i,j,1)=0;
             psm1(i,j,2)=0;
             psm1(i,j,3)=1;
        else
             psm1(i,j,1)=1;
             psm1(i,j,2)=1;
             psm1(i,j,3)=1;
        end
    end
end

%%%%%to count the pixels which diffrentiate ship from other whole area(sea)%%%%
arr=zeros(Nligfin, Ncol);
count1=0;
count2=0;
for i=1:Nligfin
    for j=1:Ncol
        if psm1(i,j,1)==0 && psm1(i,j,2)==0 && psm1(i,j,3)==0
            arr(i,j)=0;
            count1=count1+1;
        else
            arr(i,j)=1;
            count2=count2+1;
        end
    end
end

% just for checking whether we are right or not
count3=0;
count4=0;
for i=1:Nligfin
    for j=1:Ncol
        if arr(i,j)==1
            count3=count3+1;
        else
            count4=count4+1;
        end
    end
end
       
sum_ps1=0;
sum_pd1=0;
sum_pv1=0;
for i=1:Nligfin
    for j=1:Nligfin
        if arr(i,j)==1
            sum_ps1=sum_ps1+ps1(i,j);
            sum_pd1=sum_pd1+pd1(i,j);
            sum_pv1=sum_pv1+pv1(i,j);
        end
    end
end

sum_span=sum_ps1+sum_pd1+sum_pv1;
mean_ps1=sum_ps1/sum_span;
mean_pd1=sum_pd1/sum_span;
mean_pv1=sum_pv1/sum_span;
mean_span=mean_ps1+mean_pd1+mean_pv1;
% 
% 
s_ps1=0;
s_pd1=0;
s_pv1=0;
ps1=ps1./span;
pd1=(pd1)./span;
pv1=pv1./span;

for i=1:Nligfin
    for j=1:Ncol
        s_ps1=ps1(i,j)+s_ps1;
        s_pd1=pd1(i,j)+s_pd1;
        s_pv1=pv1(i,j)+s_pv1;
    end
end

m_ps1=s_ps1/(Nligfin*Ncol)
m_pd1=s_pd1/(Nligfin*Ncol)
m_pv1=s_pv1/(Nligfin*Ncol)
m_span=m_ps1+m_pd1+m_pv1













