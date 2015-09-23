
close all; 
clear all; 
clc; 

flip=0; 


%ks1, minimized 
fp=fopen('stressdistks1.stat','r'); 
NCu=35112; 
NNb=22032; 

%specify integration ranges to calculate bulk energy
r11=7;
r12=15; 
r21=27; 
r22=38; 


%{
%ks1, minimized with Liang's potential
fp=fopen('stressdistks1new.stat','r'); 
NCu=35112; 
NNb=22032; 

%specify integration ranges to calculate bulk energy
r11=7;
r12=15; 
r21=27; 
r22=38; 

%}

%{
%ks1min (i.e. with vaccancies added)
fp=fopen('stressdistksmin.stat','r'); 
NCu=35031; 
NNb=22032; 

%specify integration ranges to calculate bulk energy
r11=7;
r12=15; 
r21=27; 
r22=38; 

%}

%{
%ks1min (i.e. with vaccancies added)
fp=fopen('stressdistksminnew.stat','r'); 
NCu=35031; 
NNb=22032; 

%specify integration ranges to calculate bulk energy
r11=7;
r12=15; 
r21=27; 
r22=38; 
%}

%{
%ks2 
fp=fopen('stressdistks2.stat','r'); 
NCu=35100; 
NNb=22032; 

%specify integration ranges to calculate bulk energy
r11=52;
r12=60; 
r21=75; 
r22=83; 

flip=1; %causes us to flip data before processing b/c copper, nb layers are on wrong sides
%}

%{
%ks2new 
fp=fopen('stressdistks2new.stat','r'); 
NCu=35100; 
NNb=22032; 

%specify integration ranges to calculate bulk energy
r11=52;
r12=60; 
r21=75; 
r22=83; 

flip=1; %causes us to flip data before processing b/c copper, nb layers are on wrong sides

%}

%{
%Demk interface
fp=fopen('stressdistdemknew.stat','r'); 
NCu=30880; 
NNb=20400; 

%specify integration ranges to calculate bulk energy
r11=7;
r12=15; 
r21=27; 
r22=33;
%} 

%{
%Demk interface
fp=fopen('stressdistdemknew.stat','r'); 
NCu=30880; 
NNb=20400; 

%specify integration ranges to calculate bulk energy
r11=6;
r12=10; 
r21=27; 
r22=33; 
%}

len=fscanf(fp,'%d \n',1); 

data=textscan(fp,'%f',8*len); 
data=cell2mat(data); 

lx=fscanf(fp,'%f',1); 
ly=fscanf(fp,'%f',1); 
lz=fscanf(fp,'%f',1); 

px=zeros(1,len); 
px(1:len)=data(1:len); 
py=zeros(1,len); 
py(1:len)=data(len+1:2*len); 
pz=zeros(1,len);
pz(1:len)=data(2*len+1:3*len); 
txy=zeros(1,len); 
txy(1:len)=data(3*len+1:4*len); 
tyz=zeros(1,len); 
tyz(1:len)=data(4*len+1:5*len); 
txz=zeros(1,len); 
txz(1:len)=data(5*len+1:6*len); 
pe=zeros(1,len); 
pe(1:len)=data(6*len+1:7*len); 
counts=zeros(1,len); 
counts(1:len)=data(7*len+1:8*len); 

petemp=zeros(1,len);
countstemp=zeros(1,len); 
%flip data
if(flip==1)
    for n=1:1:len
        petemp(1,n)=pe(1,len-n+1); 
    end
    pe(1,1:len)=petemp(1,1:len); 
    for n=1:1:len
        countstemp(1,n)=counts(1,len-n+1); 
    end
    counts(1,1:len)=countstemp(1,1:len); 
end



figure; 
plot(px); 
title('PX'); 
%xlim([0,50]); 

figure; 
plot(px); 
title('PX'); 
%xlim([0,50]); 

figure; 
plot(py); 
title('PY'); 
%xlim([0,50]); 

figure; 
plot(pz); 
title('PZ'); 
%xlim([0,50]); 

figure; 
plot(txy); 
title('TXY'); 
%xlim([0,50]); 

figure; 
plot(tyz); 
title('TYZ'); 
%xlim([0,50]); 

figure; 
plot(txz); 
title('TXZ'); 
%xlim([0,50]); 

atome=zeros(1,len); 
for n=1:1:len
   if(counts(n)~=0)
        atome(n)=pe(n)/counts(n);
   end
end

figure;
plot(atome); 
title('Potential Energy'); 
%xlim([0,50]); 

xlabel('Z Distance (arbitrary)'); 
ylabel('Average Atom Potential Energy (eV)'); 

count=0; 
for n=1:1:len
    count=counts(n)+count; 
end
disp('Atom count'); 
disp(count); 

%now we want to calculate the interface energy

ecohcu=0;
count=0; 
for n=r11:1:r12
    ecohcu=ecohcu+pe(n);
    count=count+counts(n); 
end
ecohcu=ecohcu/count; 
disp('Copper cohesive energy'); 
disp(ecohcu); 

ecohnb=0;
count=0; 
for n=r21:1:r22
    ecohnb=ecohnb+pe(n);
    count=count+counts(n); 
end
ecohnb=ecohnb/count; 
disp('Nhobium cohesive energy'); 
disp(ecohnb); 

%calculate the total energy of the system and the number of atoms with each
%edge removed

%we need to start with the total number of atoms

for n=1:1:r11
    NCu=NCu-counts(n); %remove ends with free surfaces
end

s=size(counts); 
s=s(1,2); 
for n=r22:1:s
    NNb=NNb-counts(n); %remove ends with free surfaces
end

Egb=0; 
for n=(r11+1):1:(r22-1)  %sum total energy in system not including interface regions
    Egb=pe(n)+Egb; 
end

%subtract away cohesive energy
Egb=(Egb-ecohcu*NCu-ecohnb*NNb)/lx/ly; 

%Egb=Egb/(lx*10^-10)/(ly*10^-10)*1.602*10^-19*1000; %convert eV to mJ 

disp('Grain boundary energy, in eV/Ang^2'); 
disp(Egb); 

fclose(fp); 