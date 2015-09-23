
close all; 
clear all; 
clc; 

%we want to compute a representation of the elastic constant tensor Cijkl
%following a general rotation

%in addition we want to compute the forces arising from a given deformation

%this is for copper deformation
theta=0; 
phi=atan(1/2); 

T=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
%display rotation matrix   
disp(T); 

disp('T times T transpose'); 
disp(T*T'); 

disp('det T, if -1 we have reflection and rotation operation');
disp(det(T)); 

%now we need to create the elascity tensor
C=zeros(3,3,3,3); 

%for a general cubic system, there are 3 independent elastic constants,
%c11 (ciiii), c21 (ciijj), c44 (cijij), all others are zero

%elastic constants for Cu
%c11=169; 
%c12=122; 
%c44=75.3; 

%elastic constants from simulations
c11=178.6; 
c12=122.6; 
c44=80.9; 

%can also calculate for general isotropic material
%c44=(c11-c12)/2; 

for n1=1:1:3
    for n2=1:1:3
        for n3=1:1:3
            for n4=1:1:3
                if(n1==n2 && n2==n3 && n3==n4)
                    C(n1,n2,n3,n4)=c11; 
                end
                if(n1==n2 && n2~=n3 && n3==n4)
                    C(n1,n2,n3,n4)=c12; 
                end
                if(n1==n3 && n2==n4 && n3~=n2)
                    C(n1,n2,n3,n4)=c44; 
                end
                if(n1==n4 && n2==n3 && n1~=n2)
                    C(n1,n2,n3,n4)=c44; 
                end
            end
        end
    end
end

%now we want to rotate each component of the elascity tensor using the
%transformation t(i,j)

CNew=zeros(3,3,3,3); 
for n1a=1:1:3
    for n2a=1:1:3
        for n3a=1:1:3
            for n4a=1:1:3
                for n1b=1:1:3
                    for n2b=1:1:3
                        for n3b=1:1:3
                            for n4b=1:1:3
                                %two of the T operations have to be
                                %transposed
                                CNew(n1a,n2a,n3a,n4a)=T(n1b,n1a)*T(n2b,n2a)*T(n3b,n3a)*T(n4b,n4a)*C(n1b,n2b,n3b,n4b)+CNew(n1a,n2a,n3a,n4a); 
                            end
                        end
                    end
                end
            end
        end
    end
end


%{
disp1=zeros(3,3); 
force=zeros(3,3); 
disp1(1,1)=0;
disp1(1,2)=0; 
disp1(2,1)=0; 
disp1(2,2)=1; 

%calculate force using tensor operations on C
for n1a=1:1:3
    for n2a=1:1:3
        for n1b=1:1:3
            for n2b=1:1:3
                force(n1a,n2a)=C(n1a,n2a,n1b,n2b)*disp1(n1b,n2b)+force(n1a,n2a);
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
                disp2(n1a,n2a)=T(n1b,n1a)*T(n2b,n2a)*disp1(n1b,n2b)+disp2(n1a,n2a);
                force2(n1a,n2a)=T(n1b,n1a)*T(n2b,n2a)*force(n1b,n2b)+force2(n1a,n2a);
            end
        end
    end
end

%}

%set transformed matrix to new position
%C=CNew;

%after creating we want to display in the general 9 by 9 form and the 6 by
%6 form

%1=11, 4=12, 7=21 etc. 
nineform=zeros(9,9); 
nineformNew=zeros(9,9); 

for n1=1:1:3
    for n2=1:1:3
        for n3=1:1:3
            for n4=1:1:3
                
                ind1=-1; 
                ind2=-1; 
                
                if(n1==n2)
                    ind1=n1; 
                else
                    if(n1==2 && n2==3)
                        ind1=4;  
                    end
                    if(n1==3 && n2==1)
                        ind1=5;
                    end
                    if(n1==1 && n2==2)
                        ind1=6;   
                    end
                    if(n1==3 && n2==2)
                        ind1=7;   
                    end
                    if(n1==1 && n2==3)
                        ind1=8;    
                    end
                    if(n1==2 && n2==1)
                        ind1=9;    
                    end
                end
                
                if(n3==n4)
                    ind2=n3; 
                else
                    if(n3==2 && n4==3)
                        ind2=4;  
                    end
                    if(n3==3 && n4==1)
                        ind2=5;
                    end
                    if(n3==1 && n4==2)
                        ind2=6;   
                    end
                    if(n3==3 && n4==2)
                        ind2=7;   
                    end
                    if(n3==1 && n4==3)
                        ind2=8;    
                    end
                    if(n3==2 && n4==1)
                        ind2=9;    
                    end
                end
                nineformNew(ind1,ind2)=CNew(n1,n2,n3,n4); 
                nineform(ind1,ind2)=C(n1,n2,n3,n4); 
            end
        end
    end
end

disp('9 by 9 form of the elascity tensor'); 
disp(nineformNew); 

%6 by 6 form
sixformNew=zeros(6,6); 
sixform=zeros(6,6); 

for n1=1:1:3
    for n2=1:1:3
        for n3=1:1:3
            for n4=1:1:3
                
                ind1=-1; 
                ind2=-1; 
                
                if(n1==n2)
                    ind1=n1; 
                else
                    if(n1==2 && n2==3)
                        ind1=4;  
                    end
                    if(n1==3 && n2==1)
                        ind1=5;
                    end
                    if(n1==1 && n2==2)
                        ind1=6;   
                    end
                    if(n1==3 && n2==2)
                        ind1=4;   
                    end
                    if(n1==1 && n2==3)
                        ind1=5;    
                    end
                    if(n1==2 && n2==1)
                        ind1=6;    
                    end
                end
                
                if(n3==n4)
                    ind2=n3; 
                else
                    if(n3==2 && n4==3)
                        ind2=4;  
                    end
                    if(n3==3 && n4==1)
                        ind2=5;
                    end
                    if(n3==1 && n4==2)
                        ind2=6;   
                    end
                    if(n3==3 && n4==2)
                        ind2=4;   
                    end
                    if(n3==1 && n4==3)
                        ind2=5;    
                    end
                    if(n3==2 && n4==1)
                        ind2=6;    
                    end
                end
                sixform(ind1,ind2)=C(n1,n2,n3,n4); 
                sixformNew(ind1,ind2)=CNew(n1,n2,n3,n4); 
            end
        end
    end
end

disp('6 by 6 form of the elascity tensor'); 
disp(sixformNew); 




