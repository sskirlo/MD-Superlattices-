close all; 
clear all; 
clc; 

bob=fopen('cunb_gamsurf.txt','r');

inc=0.05; 
xpos=-0:inc:2; 
ypos=-0:inc:2; 
sx=size(xpos); 
sx=sx(1,2);
sy=size(ypos); 
sy=sy(1,2);
len=sx*sy; 

data=textscan(bob,'%f',len); 
data=cell2mat(data); 

%now we need to read into 2D array which we can plot surface 
z=zeros(sx,sy);
count=1; 
for n1=1:1:sx
    for n2=1:1:sy
        z(n1,n2)=data(count); 
        count=count+1; 
    end
end


%{

figure;
plot(xpoints,z(6,:));
tit=['Cross section at y=0.5 for ',a]; 
title(tit); 
xlabel('x displacement (*a0)'); 
ylabel('Grain Boundary Energy(mJ/m^2)'); 
h=gcf; 
name=[a ,'cross1', '.jpg']; 
saveas(h,name); 

figure;
plot(xpoints,z(1,:)); 
tit=['Cross section at y=0.0 for ',a]; 
title(tit);
xlabel('x displacement (*a0)'); 
ylabel('Grain Boundary Energy(mJ/m^2)'); 
h=gcf; 
name=[a ,'cross2', '.jpg']; 
saveas(h,name); 
set(gcf, 'PaperPositionMode', 'auto');
print([a,'cross2','.eps'],'-depsc2'); 

%}

figure; 
contour(xpos,ypos,z); 
tit=['Gamma Surface ']; 
title(tit);
xlabel('x displacement (*a0)'); 
ylabel('y displacement (*a0)'); 
h=gcf; 
%name=[a ,'gamma', '.jpg']; 
blah=colorbar; 
xlabel(blah,'mJ/m^2')
%saveas(h,name); 
%set(gcf, 'PaperPositionMode', 'auto');
%print([a,'gamma','.eps'],'-depsc2'); 


%identify the six absolute minimums on the surface
mins=ones(1,6); 
minx=zeros(1,6); 
miny=zeros(1,6); 


for n=1:1:6
    
    val=1000000;
    for n1=1:1:sx
        for n2=1:1:sy
            if(z(n1,n2)<val)
                mins(n)=z(n1,n2);
                minx(n)=n1;
                miny(n)=n2;
            end
        end
    end
    z(minx(n),miny(n))=0; %change to zero so next time we go around we will find different minimum
    minx(n)=xpos(minx(n));
    miny(n)=ypos(miny(n)); 
end

disp('Energy minimums and position on Gamma surface'); 
disp(mins); 
disp(minx); 
disp(miny); 

fclose(bob); 

