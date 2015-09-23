
close all; 
clear all; 
clc; 


theta=-pi/4; 
phi=0; 

TC1=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
theta=0; 
phi=-atan((2)^(1/2)); 

TC2=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
%rotate around the new <0,0,1> direction until the old <1,1,0> direction
%is in <0,1,0>, i.e. close packed direction is along y-axis

theta=pi/3; 
phi=0; 

TC3=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
TC=TC3*TC2*TC1; 

a01=3.615019; 

x1=[1,0,0]; 
x1=TC*x1.'; 
y1=[0,1,0]; 
y1=TC*y1.'; 
z1=[0,0,1]; 
z1=TC*z1.'; 

x1=x1.'; 
y1=y1.'; 
z1=z1.'; 

%this gives number of unit cells in each direction
len11=5/(2)^(1/2); 
len12=5*(3/2)^(1/2);
len13=5*(3)^(1/2);

%need to select these so that distances are bigger than simulation
%boundary, so we can be gauranteed to fill up region
%s1=ceil(len1*dot(x,[1,0,0]))+1; 
%s2=ceil(len2*dot(y,[0,1,0]))+1; 
%s3=ceil(len3*dot(z,[0,0,1]))+1; 

s11=30; 
s12=30; 
s13=30; 

%we make cell 5 times bigger than expected, and instead of modding we only 
%keep coordinates which fall within boundaries 

%define lattice vectors for each individual atom at site
%FCC lattice
llen1=4; 
vectors1=zeros(llen1,3); 
vectors1(1,:)=[0,0,0];
vectors1(2,:)=x1+y1;
vectors1(3,:)=y1+z1;
vectors1(4,:)=x1+z1;
vectors1=vectors1*1/2; 

disp(vectors1(2,:)-vectors1(3,:)); 
disp(vectors1(3,:)-vectors1(4,:)); 

coordinates1=zeros(s11*s12*s13*8*llen1,3); 
count=1;

for n1=(-s11):1:(s11-1)
    for n2=(-s12):1:(s12-1)
        for n3=(-s13):1:(s13-1)
            for n4=1:1:llen1
                pos=n1*x1+n2*y1+n3*z1+vectors1(n4,:); 
                coordinates1(count,:)=pos(:); 
                count=count+1;
            end
        end
    end
end

offset=[1,1,1]*0.01;
for n1=1:1:(s11*s12*s13*llen1*8)
    coordinates1(n1,:)=coordinates1(n1,:)+offset(:).'; 
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
keepers1=zeros(s11*s12*s13*8*llen1,3);
keepcount1=0; 
for n1=1:1:(s11*s12*s13*llen1*8)
   if(coordinates1(n1,1) >= 0 && coordinates1(n1,1) <= len11)
       if(coordinates1(n1,2) >= 0 && coordinates1(n1,2) <= len12)
           if(coordinates1(n1,3) >= 0 && coordinates1(n1,3) <= len13)
               keepcount1=keepcount1+1; 
               keepers1(keepcount1,:)=coordinates1(n1,:)*a01; 
           end
       end
   end
end

%{
figure; 
hold on; 
for n=1:1:(keepcount)
    plot3(keepers(n,1),keepers(n,2),keepers(n,3),'b*'); 
end
grid on; 
axis equal; 
axis on; 
%}

xlo=0; 
ylo=0;
zlo=0; 
xhi=len11*a01; 
yhi=len12*a01; 
zhi=len13*a01; 

fp=fopen('Cu_lattice.data','w'); 

fprintf(fp,'# Leave this line blank \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'%d atoms \n',keepcount1); 
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
for n=1:1:(keepcount1)
    fprintf(fp,'%d %d %f %f %f \n ',n,1,keepers1(n,1),keepers1(n,2),keepers1(n,3)); 
end
    
fclose(fp);
 

