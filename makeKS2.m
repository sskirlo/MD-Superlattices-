

close all; 
clear all; 
clc; 

width=67; %width of region in angstroms
gap1=2.3271; %gap between Cu-alpha and Nb layer 
gap2=2.1543; %gap between Cu-alpha and Cu layer
%gap1=30; 

anb=3.3008; 
acu=3.615019;

width1=width-mod(width,(3)^(1/2)*acu)+0.1; 
width2=width-mod(width,(2)^(1/2)*anb)+0.1; 
disp('width1'); 
disp(width1); 
disp('width2'); 
disp(width2); 

%transformation applied to Cu layer to transform it to Cu-alpha layer
Tcucua=[-(6^(1/2)-3)*(acu+anb)/acu, (3^(1/2)-2^(1/2))*(acu-anb)/acu; 
     (6-2*6^(1/2))/(2)^(1/2)+3*(2-6^(1/2))/2^(1/2)*anb/acu,(6^(1/2)-2)-(6^(1/2)-3)*anb/acu              ];


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
len11=38/(2)^(1/2); 
len12=44/2*(3/2)^(1/2);
len13=width1/a01;

%need to select these so that distances are bigger than simulation
%boundary, so we can be gauranteed to fill up region
%s1=ceil(len1*dot(x,[1,0,0]))+1; 
%s2=ceil(len2*dot(y,[0,1,0]))+1; 
%s3=ceil(len3*dot(z,[0,0,1]))+1; 

s11=60; 
s12=60; 
s13=60; 

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

%we want to apply C transformation to bottom and top of region, first need
%to filter just based on height, so won't be missing atoms in slices when
%we apply the transformations

coordinatesb=zeros(s11*s12*s13*8*llen1,3);
coordcountb=0; 
for n1=1:1:(s11*s12*s13*llen1*8)
    if(coordinates1(n1,3) >= 0 && coordinates1(n1,3) <= len13)
        coordcountb=coordcountb+1; 
        coordinatesb(coordcountb,:)=coordinates1(n1,:); 
    end
end

%now we need to select top and bottom rows and transform
%we also want to add shift to z coordinates
minr=min(coordinatesb(1:coordcountb,3)); 
maxr=max(coordinatesb(1:coordcountb,3)); 
temp=zeros(1,3); 

surfaceatoms=0; 

for n1=1:1:coordcountb
    if(coordinatesb(n1,3) >= minr-0.1 && coordinatesb(n1,3) <= minr+0.1)
        temp(1,1)=coordinatesb(n1,1)*Tcucua(1,1)+coordinatesb(n1,2)*Tcucua(1,2); 
        temp(1,2)=coordinatesb(n1,2)*Tcucua(2,2)+coordinatesb(n1,1)*Tcucua(2,1); %transform x and y coordinates
        temp(1,3)=coordinatesb(n1,3)+1/(3)^(1/2)-gap2/acu;                       %create proper gap distance
        coordinatesb(n1,:)=temp(1,:); 
        surfaceatoms=surfaceatoms+1; 
    end
    if(coordinatesb(n1,3) >= maxr-0.1 && coordinatesb(n1,3) <= maxr+0.1)
        temp(1,1)=coordinatesb(n1,1)*Tcucua(1,1)+coordinatesb(n1,2)*Tcucua(1,2); 
        temp(1,2)=coordinatesb(n1,2)*Tcucua(2,2)+coordinatesb(n1,1)*Tcucua(2,1); %transform x and y coordinates
        temp(1,3)=coordinatesb(n1,3)-1/(3)^(1/2)+gap2/acu;                       %create proper gap distance
        coordinatesb(n1,:)=temp(1,:); 
        surfaceatoms=surfaceatoms+1; 
    end
end
  
disp('Number of surface atoms transformed'); 
disp(surfaceatoms); 
  
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
keepers1=zeros(coordcountb,3);
keepcount1=0; 
for n1=1:1:coordcountb
   if(coordinatesb(n1,1) >= 0 && coordinatesb(n1,1) <= len11)
       if(coordinatesb(n1,2) >= 0 && coordinatesb(n1,2) <= len12)
           keepcount1=keepcount1+1; 
           keepers1(keepcount1,:)=coordinatesb(n1,:)*a01; 
       end
   end
end

%offset to make sure no coordinates are below 0

minz=min(keepers1(1:keepcount1,3));
for n1=1:1:keepcount1
    keepers1(n1,3)=keepers1(n1,3)-minz+0.1; 
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
len21=34*(3)^(1/2)/2; 
len22=36*(2/3)^(1/2);
len23=width2/a02;

%need to select these so that distances are bigger than simulation
%boundary, so we can be gauranteed to fill up region
%s1=ceil(len1*dot(x,[1,0,0]))+1; 
%s2=ceil(len2*dot(y,[0,1,0]))+1; 
%s3=ceil(len3*dot(z,[0,0,1]))+1; 

s21=60; 
s22=60; 
s23=60; 

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

%calculate minimum and maximum of each z component in each lattice, so we can set the gaps at the interfaces 
minCu=min(keepers1(1:keepcount1,3)); 
maxCu=max(keepers1(1:keepcount1,3)); 
minNb=min(keepers2(1:keepcount2,3)); 
maxNb=max(keepers2(1:keepcount2,3)); 

offset=[0,0,maxCu+gap1-minNb]; 
for n1=1:1:(s21*s22*s23*llen2*8)
    keepers2(n1,:)=keepers2(n1,:)+offset(:).'; 
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



xhi=97.2839; 
yhi=97.1507; 
zhi=(maxCu+gap1-minNb+maxNb+gap1-minCu); %put in gap for expected interface seperation 
xlo=0; 
ylo=0; 
zlo=0; 

disp('number of atoms'); 
disp(keepcount1+keepcount2); 

%create data file

fp=fopen('CuNb_KS22.data','w'); 

fprintf(fp,'# Leave this line blank \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'%d atoms \n',keepcount1+keepcount2); 
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
    fprintf(fp,'%d %d %f %f %f \n ',n,1,keepers1(n,1)*xhi/(len11*a01),keepers1(n,2)*yhi/(len12*a01),keepers1(n,3)); 
end
for n=1:1:(keepcount2)
    fprintf(fp,'%d %d %f %f %f \n ',n+keepcount1,2,keepers2(n,1)*xhi/(len21*a02),keepers2(n,2)*yhi/(len22*a02),keepers2(n,3)); 
end
    
fclose(fp);
 
