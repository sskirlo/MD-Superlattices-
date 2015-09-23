

close all; 
clear all; 
clc; 


a0=3.615019; 

x=[1,0,0]; 
x=x/dot(x,x)^(1/2); 
y=[0,2,-1]; 
y=y/dot(y,y)^(1/2); 
z=[0,1,2]; 
z=z/dot(z,z)^(1/2); 

%this gives number of unit cells in each direction
len1=2; 
len2=(5)^(1/2)*1; 
len3=(5)^(1/2)*1;

%need to select these so that distances are bigger than simulation
%boundary, so we can be gauranteed to fill up region
s1=ceil(len1*dot(x,[1,0,0]))+1; 
s2=ceil(len2*dot(y,[0,1,0]))+1; 
s3=ceil(len3*dot(z,[0,0,1]))+1; 

%we make cell 5 times bigger than expected, and instead of modding we only 
%keep coordinates which fall within boundaries 

%define lattice vectors for each individual atom at site
%FCC lattice
llen=4; 
vectors=zeros(llen,3); 
vectors(1,:)=[0,0,0];
vectors(2,:)=x+y;
vectors(3,:)=y+z;
vectors(4,:)=x+z;
vectors=vectors*1/2; 

coordinates=zeros(s1*s2*s3*9*llen,3); 
count=1;

for n1=0:1:(s1-1)
    for n2=0:1:(s2-1)
        for n3=0:1:(s3-1)
            for n4=1:1:llen
                pos=n1*x+n2*y+n3*z+vectors(n4,:); 
                coordinates(count,:)=pos(:); 
                count=count+1;
            end
        end
    end
end

for n1=1:1:(s1*s2*s3*llen) %copy each lattice site to 4 
    coordinates(n1+s1*s2*s3*llen,:)=coordinates(n1,:)+s2*y; 
    coordinates(n1+s1*s2*s3*2*llen,:)=coordinates(n1,:)-s2*y; 
    coordinates(n1+s1*s2*s3*3*llen,:)=coordinates(n1,:)+s3*z; 
    coordinates(n1+s1*s2*s3*4*llen,:)=coordinates(n1,:)-s3*z; 
    coordinates(n1+s1*s2*s3*5*llen,:)=coordinates(n1,:)+s2*y-s3*z; 
    coordinates(n1+s1*s2*s3*6*llen,:)=coordinates(n1,:)-s2*y-s3*z; 
    coordinates(n1+s1*s2*s3*7*llen,:)=coordinates(n1,:)+s2*y+s3*z; 
    coordinates(n1+s1*s2*s3*8*llen,:)=coordinates(n1,:)-s2*y+s3*z; 
end

offset=[1,1,1]*0.1;
for n1=1:1:(s1*s2*s3*llen*9)
    coordinates(n1,:)=coordinates(n1,:)+offset(:).'; 
end

%{
figure; 
hold on; 
for n=1:1:(s1*s2*s3*llen*9)
    plot3(coordinates(n,1),coordinates(n,2),coordinates(n,3),'b*'); 
end
grid on; 
axis equal; 
axis on; 
%}

%next we need to go through and mod each distance by the simulation box
%size
keepers=zeros(s1*s2*s3*9*llen,3);
keepcount=0; 
for n1=1:1:(s1*s2*s3*llen*9)
   if(coordinates(n1,1) >= 0 && coordinates(n1,1) <= len1)
       if(coordinates(n1,2) >= 0 && coordinates(n1,2) <= len2)
           if(coordinates(n1,3) >= 0 && coordinates(n1,3) <= len3)
               keepcount=keepcount+1; 
               keepers(keepcount,:)=coordinates(n1,:); 
           end
       end
   end
end

figure; 
hold on; 
for n=1:1:(keepcount)
    plot3(keepers(n,1),keepers(n,2),keepers(n,3),'b*'); 
end
grid on; 
axis equal; 
axis on; 


xhi=len1*a0; 
yhi=len2*a0; 
zhi=len3*a0; 
xlo=0; 
ylo=0; 
zlo=0; 

%create data file

fp=fopen('Cu_lattice.data','w'); 

fprintf(fp,'# Leave this line blank \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'%d atoms \n',keepcount); 
fprintf(fp,' \n'); 
fprintf(fp,'1 atom types \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'%f %f xlo xhi \n',xlo,xhi); 
fprintf(fp,'%f %f ylo yhi \n',ylo,yhi); 
fprintf(fp,'%f %f zlo zhi \n',zlo,zhi); 
fprintf(fp,'0.0 0.0 0.0 xy xz yz \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'Masses \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'1 63.546 \n');  
fprintf(fp,' \n'); 
fprintf(fp,'Atoms \n'); 
fprintf(fp,' \n'); 

for n=1:1:(keepcount)
    fprintf(fp,'%d %d %f %f %f \n ',n,1,keepers(n,1)*a0,keepers(n,2)*a0,keepers(n,3)*a0); 
end
    
fclose(fp); 