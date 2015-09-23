
close all; 
clear all; 
clc; 

%create Nb half of interface


theta=0; 
phi=-pi/4;  

TN1=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
% in the new <0,0,1> direction, rotate the rectangle until one of long ends
% (old <1,1,1>) is aligned with the <0,1,0> direction
   
theta=-atan((2)^(1/2)/1); 
phi=0; 

TN2=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
TN=TN2*TN1; 

a02=3.3008; 

x2=[1,0,0]; 
x2=TN*x2.'; 
y2=[0,1,0]; 
y2=TN*y2.'; 
z2=[0,0,1]; 
z2=TN*z2.'; 

x2=x2.'; 
y2=y2.'; 
z2=z2.'; 

%this gives number of unit cells in each direction
len21=8*(3)^(1/2)/2; 
len22=8*(2/3)^(1/2);
len23=8*(2)^(1/2);

%need to select these so that distances are bigger than simulation
%boundary, so we can be gauranteed to fill up region
%s1=ceil(len1*dot(x,[1,0,0]))+1; 
%s2=ceil(len2*dot(y,[0,1,0]))+1; 
%s3=ceil(len3*dot(z,[0,0,1]))+1; 

s21=30; 
s22=30; 
s23=30; 

%we make cell 5 times bigger than expected, and instead of modding we only 
%keep coordinates which fall within boundaries 

%define lattice vectors for each individual atom at site
%BCC lattice
llen2=2; 
vectors2=zeros(llen2,3); 
vectors2(1,:)=[0,0,0];
vectors2(2,:)=x2+y2+z2;
vectors2=vectors2*1/2; 


coordinates2=zeros(s21*s22*s23*8*llen2,3); 
count=1;

for n1=(-s21):1:(s21-1)
    for n2=(-s22):1:(s22-1)
        for n3=(-s23):1:(s23-1)
            for n4=1:1:llen2
                pos=n1*x2+n2*y2+n3*z2+vectors2(n4,:); 
                coordinates2(count,:)=pos(:); 
                count=count+1;
            end
        end
    end
end

%{
figure; 
hold on; 
for n=1:1:(s1*s2*s3*llen)
    plot3(coordinates(n,1),coordinates(n,2),coordinates(n,3),'b*'); 
end
grid on; 
axis equal; 
axis on; 
%}

%next we need to go through and mod each distance by the simulation box
%size
offset=[1,1,1]*0.01;
for n1=1:1:(s21*s22*s23*llen2*8)
    coordinates2(n1,:)=coordinates2(n1,:)+offset(:).'; 
end

keepers2=zeros(s21*s22*s23*8*llen2,3);
keepcount2=0; 
for n1=1:1:(s21*s22*s23*llen2*8)
   if(coordinates2(n1,1) >= 0 && coordinates2(n1,1) <= len21)
       if(coordinates2(n1,2) >= 0 && coordinates2(n1,2) <= len22)
           if(coordinates2(n1,3) >= 0 && coordinates2(n1,3) <= len23)
               keepcount2=keepcount2+1; 
               keepers2(keepcount2,:)=coordinates2(n1,:)*a02; 
           end
       end
   end
end


%{
figure; 
hold on; 
for n=1:1:(keepcount1)
    plot3(keepers1(n,1)*a01,keepers1(n,2)*a01,keepers1(n,3)*a01,'b*'); 
end
for n=1:1:(keepcount2)
    plot3(keepers2(n,1)*a02,keepers2(n,2)*a02,keepers2(n,3)*a02,'r*'); 
end
grid on; 
axis equal; 
axis on; 
%}


xlo=0; 
ylo=0;
zlo=0; 
xhi=len21*a02; 
yhi=len22*a02; 
zhi=len23*a02; 

fp=fopen('Nb_lattice.data','w'); 

fprintf(fp,'# Leave this line blank \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'%d atoms \n',keepcount2); 
fprintf(fp,' \n'); 
fprintf(fp,'3 atom types \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'%f %f xlo xhi \n',xlo,xhi); 
fprintf(fp,'%f %f ylo yhi \n',ylo,yhi); 
fprintf(fp,'%f %f zlo zhi \n',zlo,zhi); 
fprintf(fp,'0.0 0.0 0.0 xy xz yz \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'Masses \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'1 63.546 \n');  
fprintf(fp,'2 92.906 \n');
fprintf(fp,'3 4.0000 \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'Atoms \n'); 
fprintf(fp,' \n'); 

%we need to strain x and y coordinates of both lattices so they match up at
%the boundary
for n=1:1:(keepcount2)
    fprintf(fp,'%d %d %f %f %f \n ',n,2,keepers2(n,1),keepers2(n,2),keepers2(n,3)); 
end
    
fclose(fp);
 

