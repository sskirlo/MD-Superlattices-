
close all; 
clear all; 
clc; 

%specifies number of times we stack up lattice
width=2; 
label=4; %this is just a label for the file
steps=3; %put one set of steps into the interface

fp=fopen('CuNb_KS1ps.data','r');

%we want to scan in all relevant data fields
line=fgets(fp); 
disp(line); 
atomcount0=fscanf(fp,'%f ',1);
disp('atom count'); 
disp(atomcount0);
for n=1:1:4
    line=fgets(fp);
    disp(line); 
end
xlo=fscanf(fp,'%f ',1); 
xhi=fscanf(fp,'%f ',1); 
fgets(fp); 
ylo=fscanf(fp,'%f ',1); 
yhi=fscanf(fp,'%f ',1);  
fgets(fp); 
zlo=fscanf(fp,'%f ',1); 
zhi=fscanf(fp,'%f ',1); 
fgets(fp); 
disp('Simulation dimensions'); 
disp(xhi); 
disp(yhi); 
disp(zhi);
for n=1:1:9
    line=fgets(fp);
    disp(line); 
end

data=zeros(5,atomcount0); 

for n=1:1:atomcount0
    for n2=1:1:5
        data(n2,n)=fscanf(fp,'%f ',1); 
    end
end

a01=3.615019; 
a02=3.3008; 


%we need to identify the bottom two rows of the copper atoms
%and then shift additional copy until wraps back around

%search for bottom two rows of copper atoms
min1Cu=100000;
min2Cu=100000; 
cucount=0; 

for n=1:1:atomcount0
    if(data(2,n)==1) %check if its a copper atom
        if(data(5,n)<min1Cu)
            min1Cu=data(5,n); 
        end
        if(data(5,n)~=min1Cu && data(5,n)<min2Cu)
            min2Cu=data(5,n);
        end
        cucount=cucount+1; 
    end
end

min1Nb=100000;
min2Nb=100000; 
nbcount=0; 

for n=1:1:atomcount0
    if(data(2,n)==2) %check if its a copper atom
        if(data(5,n)<min1Nb)
            min1Nb=data(5,n); 
        end
        if(data(5,n)~=min1Nb && data(5,n)<min2Nb)
            min2Nb=data(5,n);
        end
        nbcount=nbcount+1; 
    end
end

max1Cu=-100000;
max2Cu=-100000; 
cucount=0; 

for n=1:1:atomcount0
    if(data(2,n)==1) %check if its a copper atom
        if(data(5,n)>max1Cu)
            max1Cu=data(5,n); 
        end
        if(data(5,n)~=max1Cu && data(5,n)>max2Cu)
            max2Cu=data(5,n);
        end
    end
end

max1Nb=-100000;
max2Nb=-100000; 
cucount=0; 

for n=1:1:atomcount0
    if(data(2,n)==2) %check if its a copper atom
        if(data(5,n)>max1Nb)
            max1Nb=data(5,n); 
        end
        if(data(5,n)~=max1Nb && data(5,n)>max2Nb)
            max2Nb=data(5,n);
        end
    end
end

%this is step where we duplicate lattice
%final duplication step occurs when creating step

data2=data; 
datatemp=zeros(5,atomcount0); 
for n2=1:1:(width)
    for n=1:1:atomcount0
        datatemp(:,n)=data(:,n); 
        datatemp(1,n)=datatemp(1,n)+atomcount0*n2;  %displace number
        datatemp(4,n)=data(4,n)+yhi*n2;             %displace dimension
    end
    data2=[data2,datatemp]; %we keep on adding new atom coordinates 
end

data=data2; 
atomcount=atomcount0*width; 

%we need to remove sections with min and max for Cu and Nb layers

%count number of nb and cu atoms in top and bottom layers
countbotcu=0; 
counttopcu=0;
countbotnb=0; 
counttopnb=0;
for n=1:1:atomcount
    if(data(5,n)<(min2Cu+0.1) && data(5,n)>(min1Cu-0.1) )
        countbotcu=countbotcu+1; 
    end
    if(data(5,n)<(min2Nb+0.1) && data(5,n)>(min1Nb-0.1) )
        countbotnb=countbotnb+1; 
    end
    if(data(5,n)<(max1Cu+0.1) && data(5,n)>(max2Cu-0.1) )
        counttopcu=counttopcu+1; 
    end
    if(data(5,n)<(max1Nb+0.1) && data(5,n)>(max2Nb-0.1) )
        counttopnb=counttopnb+1; 
    end
end

botCu=zeros(5,countbotcu); 
topCu=zeros(5,counttopcu); 
botNb=zeros(5,countbotnb); 
topNb=zeros(5,counttopnb); 

%we create lists of atoms in bottom and top layers which we can process
%more rapidly than all of the atoms in the system

countbotcu=0; 
counttopcu=0;
countbotnb=0; 
counttopnb=0;
for n=1:1:atomcount
    if(data(5,n)<(min2Cu+0.1) && data(5,n)>(min1Cu-0.1) )
        countbotcu=countbotcu+1; 
        botCu(:,countbotcu)=data(:,n); 
        data(2,n)=-1; %label as removed atom
    end
    if(data(5,n)<(min2Nb+0.1) && data(5,n)>(min1Nb-0.1) )
        countbotnb=countbotnb+1; 
        botNb(:,countbotnb)=data(:,n); 
        data(2,n)=-1; %label as removed atom
    end
    if(data(5,n)<(max1Cu+0.1) && data(5,n)>(max2Cu-0.1) )
        counttopcu=counttopcu+1; 
        topCu(:,counttopcu)=data(:,n); 
        data(2,n)=-1; %label as removed atom
    end
    if(data(5,n)<(max1Nb+0.1) && data(5,n)>(max2Nb-0.1) )
        counttopnb=counttopnb+1; 
        topNb(:,counttopnb)=data(:,n); 
        data(2,n)=-1; %label as removed atom
    end
end

disp('min and max coordinates'); 
disp(min1Cu); 
disp(min2Cu); 
disp(max1Cu); 
disp(max2Cu); 
disp(min1Nb); 
disp(min2Nb); 
disp(max1Nb); 
disp(max2Nb); 

disp('total atoms subtracted from main lists'); 
disp(countbotcu+countbotnb+counttopcu+counttopnb); 

cs=a01*(6)^(1/2); %we want to make sure we get an complete unit cells of Ni on each terrace
csad=0.02; 

countcurm=0; 
%now we need to remove atoms from these restricted lists
%remove atoms from top Cu list
for n=1:1:counttopcu
    for n2=1:1:steps
        if(topCu(4,n)>(ceil(yhi*width/(steps*2)*(n2*2-1)/cs)*cs+csad) && topCu(4,n)<(ceil(yhi*width/(steps*2)*(n2*2)/cs)*cs+csad)) %make sure in the valid step region
            if(topCu(5,n)>(max2Cu+0.1) && topCu(5,n)< (max1Cu+0.1))
                topCu(2,n)=-1; %label as removed atom
                countcurm=countcurm+1; 
            end
        end
    end
end

countnbrm=0; 
%now we need to remove atoms from these restricted lists
%remove atoms from top Nb list
for n=1:1:counttopnb
    for n2=1:1:steps
        if(topNb(4,n)>(ceil(yhi*width/(steps*2)*(n2*2-1)/cs)*cs+csad) && topNb(4,n)<(ceil(yhi*width/(steps*2)*(n2*2)/cs)*cs+csad)) %make sure in the valid step region
            if(topNb(5,n)>(max2Nb+0.1) && topNb(5,n)< (max1Nb+0.1))
                topNb(2,n)=-1; %label as removed atom
                countnbrm=countnbrm+1; 
            end
        end
    end
end

disp('Removed atoms'); 
disp('Cu removed');
disp(countcurm); 
disp('Nb removed'); 
disp(countnbrm); 


%define cu nearest neighbor vectors
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



x1=[1,0,0]; 
x1=TC*x1.'; 
y1=[0,1,0]; 
y1=TC*y1.'; 
z1=[0,0,1]; 
z1=TC*z1.'; 

x1=x1.'*a01; 
y1=y1.'*a01; 
z1=z1.'*a01; 

llen1=6; 
vectors1=zeros(llen1,3); 
vectors1(1,:)=x1+y1;
vectors1(2,:)=y1+z1;
vectors1(3,:)=x1+z1;
vectors1=vectors1*1/2; 

%we need to add these six vectors to bottom copper list and see which ones
%end up in region where Nb atoms where removed

addcus=zeros(5,countbotcu*4); 
addcucount=1; 



for n=1:1:countbotcu
    for n3=1:1:3
        addcus(2,n)=1; %label as Cu atom
        addcus(3:5,addcucount)=botCu(3:5,n)-vectors1(n3,1:3).'; 
        for n2=1:1:steps
            if(addcus(4,addcucount)>(ceil(yhi*width/(steps*2)*(n2*2-1)/cs)*cs+csad) && addcus(4,addcucount)<(ceil(yhi*width/(steps*2)*(n2*2)/cs)*cs+csad)) %make sure in the valid step region
                if(addcus(5,addcucount)>(max2Nb+0.5-zhi) && addcus(5,addcucount)< (max1Nb+0.5-zhi)) %wrap around to bottom
                    %we need to wrap around the boundaries
                    if(addcus(3,addcucount)>xhi)
                        addcus(3,addcucount)=addcus(3,addcucount)-xhi; 
                    end
                    if(addcus(3,addcucount)<0)
                        addcus(3,addcucount)=addcus(3,addcucount)+xhi; 
                    end
                    addcucount=addcucount+1; 
                end
            end
        end
    end
end
addcucount=addcucount-1; 

%we need to wrap around coordinates, and eliminate redundant coordinates

finalcucount=0;
finalcu=zeros(5,countbotcu); 

for n=1:1:addcucount
    val=-1; 
    %need to check if atom coordinate is already in list
    for n2=1:1:finalcucount
        if((finalcu(3,n2)>addcus(3,n)-0.01) && (finalcu(3,n2)<addcus(3,n)+0.01) && (finalcu(4,n2)>addcus(4,n)-0.01) && (finalcu(4,n2)<addcus(4,n)+0.01))
            val=1; 
        end
    end
    if(val==-1)
        finalcucount=finalcucount+1; 
        finalcu(:,finalcucount)=addcus(:,n); 
    end
end

disp('add Cu count (initially added)'); 
disp(addcucount); 

disp('final Cu count (after redundant atoms eliminated)'); 
disp(finalcucount); 


%Nb atom count

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


x2=[1,0,0]; 
x2=TN*x2.'; 
y2=[0,1,0]; 
y2=TN*y2.'; 
z2=[0,0,1]; 
z2=TN*z2.'; 

x2=x2.'*a02; 
y2=y2.'*a02; 
z2=z2.'*a02;

llen2=4; 
vectors2=zeros(llen2,3); 
vectors2(1,:)=x2+y2+z2;
vectors2(2,:)=-x2+y2+z2;
vectors2(3,:)=x2-y2-z2;
vectors2(4,:)=-x2-y2-z2;
vectors2=vectors2*1/2; 

addnbs=zeros(5,countbotnb*4); 
addnbcount=1; 



for n=1:1:countbotnb
    for n3=1:1:4
        addnbs(2,n)=2; %label as Nb atom
        addnbs(3:5,addnbcount)=botNb(3:5,n)+vectors2(n3,1:3).'; 
        for n2=1:1:steps
            if(addnbs(4,addnbcount)>(ceil(yhi*width/(steps*2)*(n2*2-1)/cs)*cs+csad) && addnbs(4,addnbcount)<(ceil(yhi*width/(steps*2)*(n2*2)/cs)*cs+csad)) %make sure in the valid step region
                if(addnbs(5,addnbcount)>(max2Cu+0.5) && addnbs(5,addnbcount)< (max1Cu+0.5)) %wrap around to bottom
                    %we need to wrap around the boundaries
                    if(addnbs(3,addnbcount)>xhi)
                        addnbs(3,addnbcount)=addnbs(3,addnbcount)-xhi; 
                    end
                    if(addnbs(3,addnbcount)<0)
                        addnbs(3,addnbcount)=addnbs(3,addnbcount)+xhi; 
                    end
                    addnbcount=addnbcount+1; 
                end
            end
        end
    end
end
addnbcount=addnbcount-1; 

%we need to wrap around coordinates, and eliminate redundant coordinates

finalnbcount=0;
finalnb=zeros(5,countbotnb*2); 

for n=1:1:addnbcount
    val=-1; 
    %need to check if atom coordinate is already in list
    for n2=1:1:finalnbcount
        if((finalnb(3,n2)>addnbs(3,n)-0.01) && (finalnb(3,n2)<addnbs(3,n)+0.01) && (finalnb(4,n2)>addnbs(4,n)-0.01) && (finalnb(4,n2)<addnbs(4,n)+0.01))
            val=1; 
        end
    end
    if(val==-1)
        finalnbcount=finalnbcount+1; 
        finalnb(:,finalnbcount)=addnbs(:,n); 
    end
end

disp('add Nb count (initially)'); 
disp(addnbcount); 

disp('final Nb count (after redundant atoms eliminated)'); 
disp(finalnbcount); 



%we sort into Cu and Nb lists

dataCu=zeros(5,cucount*width); 
dataNb=zeros(5,nbcount*width);
count1=0; 
count2=0; 

mins=ones(1,3)*10000; 
maxs=ones(1,3)*-10000; 

for n=1:1:atomcount
    
    for n2=1:1:3
        if(mins(1,n2)>data(n2+2,n))
            mins(1,n2)=data(n2+2,n); 
        end
        if(maxs(1,n2)<data(n2+2,n))
            maxs(1,n2)=data(n2+2,n); 
        end
    end
        
    if(data(2,n)==1)
        count1=count1+1; 
        dataCu(:,count1)=data(:,n); 
    end
    if(data(2,n)==2)
        count2=count2+1; 
        dataNb(:,count2)=data(:,n); 
    end
end

for n=1:1:counttopnb
    if(topNb(2,n)==1)
        count1=count1+1; 
        dataCu(:,count1)=topNb(:,n); 
    end
    if(topNb(2,n)==2)
        count2=count2+1; 
        dataNb(:,count2)=topNb(:,n); 
    end
end
for n=1:1:counttopcu
    if(topCu(2,n)==1)
        count1=count1+1; 
        dataCu(:,count1)=topCu(:,n); 
    end
    if(topCu(2,n)==2)
        count2=count2+1; 
        dataNb(:,count2)=topCu(:,n); 
    end
end
for n=1:1:countbotnb
    if(botNb(2,n)==1)
        count1=count1+1; 
        dataCu(:,count1)=botNb(:,n); 
    end
    if(botNb(2,n)==2)
        count2=count2+1; 
        dataNb(:,count2)=botNb(:,n); 
    end
end
for n=1:1:countbotcu
    if(botCu(2,n)==1)
        count1=count1+1; 
        dataCu(:,count1)=botCu(:,n); 
    end
    if(botCu(2,n)==2)
        count2=count2+1; 
        dataNb(:,count2)=botCu(:,n); 
    end
end
%add in additional rows
for n=1:1:finalnbcount
    if(finalnb(2,n)==2)
        count2=count2+1; 
        dataNb(:,count2)=finalnb(:,n); 
    end
end
for n=1:1:finalcucount
    if(finalcu(2,n)==1)
        count1=count1+1; 
        dataCu(:,count1)=finalcu(:,n); 
    end
end

disp('Max dimensions and min dimensions');
disp(maxs); 
disp(mins); 

disp('Compare counts'); 
disp(atomcount); 
disp(count1+count2); 

a=sprintf('%d',label); 
fp=fopen(['KS1steps',a,'.data'],'w'); 

fprintf(fp,'# Leave this line blank \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'%d atoms \n',count1+count2); 
fprintf(fp,' \n'); 
fprintf(fp,'3 atom types \n'); 
fprintf(fp,' \n'); 
fprintf(fp,'%f %f xlo xhi \n',xlo,xhi); 
fprintf(fp,'%f %f ylo yhi \n',ylo,yhi*width); 
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

idcount=1; 
for n=1:1:count1
    fprintf(fp,'%d %d %f %f %f \n',idcount,dataCu(2,n),dataCu(3,n),dataCu(4,n),dataCu(5,n)); 
    idcount=idcount+1; 
end
for n=1:1:count2
    fprintf(fp,'%d %d %f %f %f \n',idcount,dataNb(2,n),dataNb(3,n),dataNb(4,n),dataNb(5,n)); 
    idcount=idcount+1; 
end

fclose(fp); 

%{
figure; 
hold on; 
for n=1:10:cucount*2^(width)
    plot3(dataCu(3,n),dataCu(4,n),dataCu(5,n),'b*'); 
end
for n=1:10:nbcount*2^(width)
    plot3(dataNb(3,n),dataNb(4,n),dataNb(5,n),'r*'); 
end
%}

%xlim([mins(1,1),maxs(1,1)]); 
%ylim([mins(1,2),maxs(1,2)]); 
%zlim([mins(1,3),maxs(1,3)]); 



