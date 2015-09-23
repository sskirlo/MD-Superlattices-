
close all; 
clear all; 
clc; 

fp=fopen('trial3.txt','r'); 

readin=' '; 
for n=1:1:15
    readin=[readin,' %f']; 
end
data=textscan(fp,readin); 
data=cell2mat(data); 

sl=size(data); 
sl=sl(1,1); 

%data columns step pe xhi yhi zhi xy yz xz pxx pyy pzz pxy pyz pxz

%lets calculate elastic constants

%use these as short hand for looking up components
xhi=3;
yhi=4;
zhi=5; 
xy=6; 
yz=7; 
xz=8; 
pxx=9; 
pyy=10; 
pzz=11; 
pxy=12; 
pyz=13; 
pxz=14; 

%we want to calculate the effective elastic constant c1111 for the rotated
%system need to project strain in system into this (hopefully will just yield strain in this direction)

%following this we need to project the pressure tensor into this direction
%the slope of the project pressure vs the strain will give c1111. 

%need negative of original angles to rotate deformation strain back into 1
%direction

%trial 1
%theta=0; 
%phi=0; 

%trial 2
%theta=pi/4; 
%phi=0; 

%trial 3
theta=pi/3; 
phi=pi/3; 

T=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   




strain=zeros(6,sl); 
stress=zeros(6,sl);
%extract first three components of resultant fields
for n2=1:1:3
    for n=1:1:sl
        strain(n2,n)=(data(n,xhi+n2-1)-data(1,xhi+n2-1))/data(1,xhi+n2-1); 
        stress(n2,n)=data(n,pxx+n2-1)/10; 
    end
end
%need to last 6 differently b/c shear components require different
%normalization
for n2=4:1:6
    for n=1:1:sl
        if(n2==4)
            len=data(1,zhi);
        end
        if(n2==5)
            len=data(1,xhi);
        end
        if(n2==6)
            len=data(1,yhi);
        end
        strain(n2,n)=(data(n,xhi+n2-1))/len*1/2; %need to divide by 2 also b/c shear strain
        stress(n2,n)=data(n,pxx+n2-1)/10; 
    end
end

%now we want to calculate projections of stress and strain onto 1,1
%direction

%need to read strain, stress into tensors so can transform properly

strainf=zeros(1,sl); 
stressf=zeros(1,sl); 

for n=1:1:sl

    displ=zeros(3,3); 
    force=zeros(3,3); 
    for n1=1:1:3
        for n2=1:1:3
            if(n1==n2)
                displ(n1,n1)=strain(n1,n); 
                force(n1,n1)=stress(n1,n); 
            else
                if((n1==1 && n2==2) || (n1==2 && n2==1))
                    displ(n1,n2)=strain(4,n); 
                    force(n1,n2)=stress(4,n); 
                end
                if((n1==2 && n2==3) || (n1==3 && n2==2))
                    displ(n1,n2)=strain(5,n); 
                    force(n1,n2)=stress(5,n); 
                end
                if((n1==1 && n2==3) || (n1==3 && n2==1))
                    displ(n1,n2)=strain(6,n); 
                    force(n1,n2)=stress(6,n); 
                end
            end
        end
    end
    disp2=zeros(3,3); 
    force2=zeros(3,3); 
    for n1a=1:1:3
        for n2a=1:1:3
            for n1b=1:1:3
                for n2b=1:1:3
                    disp2(n1a,n2a)=T(n1b,n1a)*T(n2b,n2a)*displ(n1b,n2b)+disp2(n1a,n2a);
                    force2(n1a,n2a)=T(n1b,n1a)*T(n2b,n2a)*force(n1b,n2b)+force2(n1a,n2a);
                end
            end
        end
    end
    
    strainf(n)=disp2(2,2); 
    stressf(n)=force2(2,2); 
    %now that we have each tensor, we need to transform the components 


end

c1=polyfit(strainf,stressf,1); 
x1=-0.001:0.0001:0.001; 
s=size(x1); 
s=s(1,2);
y1=zeros(1,s);
for n=1:1:s
    y1(n)=c1(1,1)*x1(n)+c1(1,2); 
end

figure; 
plot(strainf,stressf,'b'); 
hold on; 
plot(x1,y1,'r'); 
title('Pxx vs. x strain');

%disp('c1111 in rotated system');  
disp(c1(1,1)); 

fclose(fp); 