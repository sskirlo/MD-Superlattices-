

close all; 
clear all; 
clc; 


width=40; %width of region in angstroms
gap=2.4203;
%gap=5; 
steps=3; 
scale=2; 

stepz=zeros(1,steps*2); 

%we want to create the periodic step
for n=1:2:(steps*2)
    stepz(1,n)=1; 
end
    
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
len12=44/2*(3/2)^(1/2)*scale;
len13=width/a01;

%need to select these so that distances are bigger than simulation
%boundary, so we can be gauranteed to fill up region
%s1=ceil(len1*dot(x,[1,0,0]))+1; 
%s2=ceil(len2*dot(y,[0,1,0]))+1; 
%s3=ceil(len3*dot(z,[0,0,1]))+1; 

s11=60*scale/2; 
s12=60*scale/2; 
s13=60*scale/2; 

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
maxCu=-10000; 
minCu=10000; 
for n1=1:1:(s11*s12*s13*llen1*8)
   if(coordinates1(n1,1) >= 0 && coordinates1(n1,1) <= len11)
       true=0; 
       for n2=1:1:(steps*2)  %depending on y coordinate, we need to apply different restrictions to x coordinate
           if(coordinates1(n1,2)>=len12*(n2-1)/(steps*2) && coordinates1(n1,2) <=len12*(n2)/(steps*2) )
               if(coordinates1(n1,3) >= 0+stepz(1,n2)*0.5 && coordinates1(n1,3) <= len13+stepz(1,n2)*0.5)
                   true=1; 
                   %need to determine min and max in a given slice
                   if(n2==1)
                       if(coordinates1(n1,3)>maxCu)
                           maxCu=coordinates1(n1,3); 
                       end
                       if(coordinates1(n1,3)<minCu)
                           minCu=coordinates1(n1,3); 
                       end
                   end
                   
               end
           end
       end
       if(true==1)
           keepcount1=keepcount1+1; 
           keepers1(keepcount1,:)=coordinates1(n1,:)*a01; 
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
len22=36*(2/3)^(1/2)*scale;
len23=width/a02;

%need to select these so that distances are bigger than simulation
%boundary, so we can be gauranteed to fill up region
%s1=ceil(len1*dot(x,[1,0,0]))+1; 
%s2=ceil(len2*dot(y,[0,1,0]))+1; 
%s3=ceil(len3*dot(z,[0,0,1]))+1; 

s21=80*scale/2; 
s22=80*scale/2; 
s23=80*scale/2; 

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
minNb=10000; 
maxNb=-10000; 
for n1=1:1:(s21*s22*s23*llen2*8)
   if(coordinates2(n1,1) >= 0 && coordinates2(n1,1) <= len21)
       true=0; 
       for n2=1:1:(steps*2)  %depending on y coordinate, we need to apply different restrictions to x coordinate
           if(coordinates2(n1,2)>=len22*(n2-1)/(steps*2) && coordinates2(n1,2) <=len22*(n2)/(steps*2) )
               if(coordinates2(n1,3) >= 0+stepz(1,n2)*0.7 && coordinates2(n1,3) <= len23+stepz(1,n2)*0.7)
                   true=1; 
                   %need to determine min and max in a given slice
                   if(n2==1)
                       if(coordinates2(n1,3)>maxNb)
                           maxNb=coordinates2(n1,3); 
                       end
                       if(coordinates2(n1,3)<minNb)
                           minNb=coordinates2(n1,3); 
                       end
                   end
               end
           end
       end
       if(true==1)
           keepcount2=keepcount2+1; 
           keepers2(keepcount2,:)=coordinates2(n1,:)*a02; 
       end
   end
end

%calculate minimum and maximum of each z component in each lattice, so we can set the gaps at the interfaces 
%minCu=min(keepers1(1:keepcount1,3)); 
%maxCu=max(keepers1(1:keepcount1,3)); 
%minNb=min(keepers2(1:keepcount2,3)); 
%maxNb=max(keepers2(1:keepcount2,3)); 

minCu=minCu*a01; 
maxCu=maxCu*a01; 
minNb=minNb*a02; 
maxNb=maxNb*a02; 

offset=[0,0,maxCu+gap-minNb]; 
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
yhi=97.1507*scale; 
zhi=(maxCu+gap-minNb+maxNb+gap-minCu); %put in gap for expected interface seperation 
xlo=0; 
ylo=0; 
zlo=0; 

%create data file

a=sprintf('%d',steps); 
fp=fopen(['KS1steps',a,'.data'],'w'); 

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
 
