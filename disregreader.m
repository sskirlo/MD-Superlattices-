
close all; 
clear all; 
clc; 

fp=fopen('disregistry1.txt','r'); 

data=textscan(fp,'%f %f %f %f %f %f %f %f');
data=cell2mat(data); 

s=size(data); 
s=s(1,1); 

%we want to determine the x and y coordinates of each vector along with the
%2D position to plot

% id1 1 id2 2 vec1 3 vec2 4 vec3 5 pos1 6 pos2 7 pos3 8

x=data(:,6); 
y=data(:,7); 
u=data(:,3); 
v=data(:,4); 

figure; 
quiver(x,y,u,v);
title('Interface Disregistry with Cu_{\alpha} reference structure'); 
xlabel('x dimension (A)'); 
xlabel('y dimension (A)'); 
axis equal; 

mag=zeros(1,s); 
for n=1:1:s
    mag(n)=(u(n)^2+v(n)^2)^(1/2); 
end

scale=2; 
%we want to hash in magnitudes of the vectors
len1=40; 
len2=40; 
magsquare=zeros(len1,len2); 
ids=ones(len1,len2)*-1; 

filled=0; 
for n=1:1:s
    x1=ceil(x(n)/scale);
    y1=ceil(y(n)/scale);
    if(x1<(len1+1) && x1>0 && y1>0 && y1<(len2+1))
        magsquare(x1,y1)=mag(n); 
        ids(x1,y1)=n;
        filled=filled+1; 
    end
end

disp('filled spots of len1*len2 '); 
disp(filled); 


%now we need to go through and fill in squares with zero value
%(or with id -1)
%we do this by looking at neighbors and picking ones closest in proximity

%{
for n2=1:1:len2
    for n1=1:1:len1
        if(ids(n1,n2)==-1)
            %find neighbors
            %calculate cell center
            x0=n1*scale-0.5; 
            y0=n2*scale-0.5;
            distance=1000; 
            keeper=-1; 
            for n3=-1:1:1
                for n4=-1:1:1
                    if((n1+n3)>0 && (n1+n3)<(len1+1) && (n2+n4)>0 && (n2+n4)<(len2+1) )
                        %check distance, if there
                        if(ids(n1+n3,n2+n4)~=-1)
                            id=ids(n1+n3,n2+n4); 
                            if( ((x(id)-x0)^2+(y(id)-y0)^2)^(1/2) < distance)
                                distance=((x(id)-x0)^2+(y(id)-y0)^2)^(1/2); 
                                keeper=id; 
                            end
                        end
                    end
                end
            end
            if(keeper~=-1)
                magsquare(n1,n2)=mag(keeper); 
                ids(n1,n2)=keeper; 
            end
        end
    end
end


figure; 
contour(magsquare); 
[vx,vy]=gradient(magsquare); 

%now we want to find the direction of steepest decent, just average
%gradient over window in middle, this gives the direction

sumx=0;
sumy=0; 
counts=0; 
for n1=1:1:len1
    for n2=1:1:len2
        sumx=sumx+vx(n1,n2); 
        sumy=sumy+vy(n1,n2); 
        counts=counts+1; 
    end
end

disp('normal vector to dislocation direction '); 
disp(sumx/counts); 
disp(sumy/counts); 

if(sumx<0)
    sumx=sumx*-1; 
    sumy=sumy*-1; 
end

%}

%eigenvectors from analytic results
 %-0.5000    0.4611
 %0.8660   -0.8874

normx=0.8874; 
normy=0.4611;

norm=[normx,normy]; 
norm=norm/dot(norm,norm)^(1/2); 

xs=zeros(1,80); 
ys=zeros(1,80); 
for n=1:1:80
    xs(n)=norm(1,1)*n; 
    ys(n)=norm(1,2)*n; 
end
    
hold on; 
plot(xs,ys,'r'); 
axis equal; 

%now we can actually try to visualize the dislocation cores

par=[normy,-normx]; 

screw=zeros(1,s); 
edge=zeros(1,s);
position=zeros(1,s); %this gives the distance along 

for n=1:1:s
    screw(n)=dot(par,[u(n),v(n)]); 
    edge(n)=dot(norm,[u(n),v(n)]); 
    position(n)=dot(norm,[x(n),y(n)]); 
end

figure; 
plot(position,screw,'r.');
hold on; 
plot(position,edge,'b.'); 
title('Edge and screw components of dislocation'); 
xlabel('Perpindicular distance to dislocations (A)');
ylabel('Disregistry component (A)'); 

