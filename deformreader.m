
close all; 
clear all; 
clc; 

fp=fopen('cuc11.txt','r'); 

readin=' '; 
for n=1:1:18
    readin=[readin,' %f']; 
end
data=textscan(fp,readin); 
data=cell2mat(data); 

s=size(data); 
s=s(1,1); 

%data columns step pe xlo xhi ylo yhi zlo zhi xy yz xz press pxx pyy pzz pxy pyz pxz

%lets calculate elastic constants

%use these as short hand for looking up components
xhi=4;
yhi=6;
zhi=8; 
xy=9; 
yz=10; 
xz=11; 
pxx=13; 
pyy=14; 
pzz=15; 
pxy=16; 
pyz=17; 
pxz=18; 


start=data(1,xhi); 

strain=zeros(1,s); 
stress=zeros(1,s); 
for n=1:1:s
    strain(n)=(data(n,xhi)-start)/start; 
    stress1(n)=data(n,pxx)/10; 
    stress2(n)=data(n,pyy)/10; 
end

c1=polyfit(strain,stress1,1); 
x1=-0.001:0.0001:0.001; 
s=size(x1); 
s=s(1,2);
y1=zeros(1,s);
for n=1:1:s
    y1(n)=c1(1,1)*x1(n)+c1(1,2); 
end

figure; 
plot(strain,stress1,'b'); 
hold on; 
plot(x1,y1,'r'); 
title('Pxx vs. x strain')

c2=polyfit(strain,stress2,1); 
x2=-0.001:0.0001:0.001; 
s=size(x2); 
s=s(1,2);
y2=zeros(1,s);
for n=1:1:s
    y2(n)=c2(1,1)*x2(n)+c2(1,2); 
end

figure; 
plot(strain,stress2,'b'); 
hold on; 
plot(x2,y2,'r'); 
title('Pyy vs. x strain'); 

fclose(fp); 

%open another file to get elastic constant from shear test

fp=fopen('cuc44.txt','r'); 

data=textscan(fp,readin); 
data=cell2mat(data); 

s=size(data); 
s=s(1,1); 

strain2=zeros(1,s); 
stress3=zeros(1,s); 
for n=1:1:s
    strain2(n)=(data(n,xy)-data(1,xy))/data(1,zhi);  %we don't have symmetric component  
    stress3(n)=data(n,pxy)/10; 
end

c3=polyfit(strain2,stress3,1); 
x3=-0.0005:0.00001:0.0005; 
s=size(x3); 
s=s(1,2);
y3=zeros(1,s);
for n=1:1:s
    y3(n)=c3(1,1)*x3(n)+c3(1,2); 
end

figure; 
plot(strain2,stress3,'b'); 
hold on; 
plot(x3,y3,'r'); 
title('Pxy vs. xy strain'); 


disp('c11 in MPa'); 
disp(c1(1,1)); 
disp('c12 in MPa'); 
disp(c2(1,1)); 
disp('c44 in MPa'); 
disp(c3(1,1)); 



fclose(fp); 