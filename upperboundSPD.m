
close all; 
clear all; 
clc; 

%we want to calculate the upper bound for the elascity tensor of the Cu Nb
%bilayer 

%first we need to transform each elascity tensor into the orientation that
%the two crystals are joined in


%Cu transformations

theta=-pi/4; 
phi=0; 

%rotate <-1,1,0> direction to y axis
TC1=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
theta=0; 
phi=-atan(1/(2)^(1/2)); 

%rotate <1,1,-1> direction up to z plane
TC2=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
  
   
%Nb transformations   
   
theta=-pi/4; 
phi=0; 

%rotate <-1,1,0> direction to y axis
TN1=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
theta=0; 
phi=-atan(1/(2)^(1/2)); 

%rotate <1,1,-1> direction up to z plane
TN2=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
theta=pi/2; 
phi=0; 

%do a final 90 degree rotation so we can get <1,-1,0> direction along x
%axis and <1,1,-1> direction along y-axis
TN3=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ];   

   
%we want to check certain features of each set of transformations to make
%sure they are preforming as expected

%we want to calculate a complete 3 by 3 matrix which does all 3
%transformations
TC=TC1*TC2; 
TN=TN1*TN2*TN3;  

%checks for Cu
vec1c=[-1,1,-2]/6^(1/2); 
vec2c=[1,1,0]/2^(1/2); 
vec3c=[-1,1,1]/3^(1/2); 

res1=TC.'*vec1c.'; 
res2=TC.'*vec2c.'; 
res3=TC.'*vec3c.';  

disp('vec 1 <-1,1,-2>'); 
disp(res1); 
disp('vec 2 <1,-1,0>'); 
disp(res2); 
disp('vec 3 <1,1,-1>'); 
disp(res3); 


%checks for Nb
vec1n=[-1,1,-2]/6^(1/2); 
vec2n=[1,1,0]/2^(1/2); 
vec3n=[-1,1,1]/3^(1/2); 

res1=TN.'*vec1n.'; 
res2=TN.'*vec2n.'; 
res3=TN.'*vec3n.'; 

disp('vec 1 <1,1,-2>'); 
disp(res1); 
disp('vec 2 <1,-1,0>'); 
disp(res2); 
disp('vec 3 <1,1,-1>'); 
disp(res3); 

%now we need to define both the Cu and Nb stress tensors and preform the
%rotations on them necessary to represent the two crystals at the interface

CCu=zeros(3,3,3,3); 

%elastic constants from simulations for Cu
%c11=178.6; 
%c12=122.6; 
%c44=80.9; 
c11=179; 
c12=123; 
c44=81; 

for n1=1:1:3
    for n2=1:1:3
        for n3=1:1:3
            for n4=1:1:3
                if(n1==n2 && n2==n3 && n3==n4)
                    CCu(n1,n2,n3,n4)=c11; 
                end
                if(n1==n2 && n2~=n3 && n3==n4)
                    CCu(n1,n2,n3,n4)=c12; 
                end
                if(n1==n3 && n2==n4 && n3~=n2)
                    CCu(n1,n2,n3,n4)=c44; 
                end
                if(n1==n4 && n2==n3 && n1~=n2)
                    CCu(n1,n2,n3,n4)=c44; 
                end
            end
        end
    end
end

CNb=zeros(3,3,3,3); 

%elastic constants from simulations for Cu
c11=245.8; 
c12=133.8; 
c44=28.8; 

for n1=1:1:3
    for n2=1:1:3
        for n3=1:1:3
            for n4=1:1:3
                if(n1==n2 && n2==n3 && n3==n4)
                    CNb(n1,n2,n3,n4)=c11; 
                end
                if(n1==n2 && n2~=n3 && n3==n4)
                    CNb(n1,n2,n3,n4)=c12; 
                end
                if(n1==n3 && n2==n4 && n3~=n2)
                    CNb(n1,n2,n3,n4)=c44; 
                end
                if(n1==n4 && n2==n3 && n1~=n2)
                    CNb(n1,n2,n3,n4)=c44; 
                end
            end
        end
    end
end

CNCu=zeros(3,3,3,3); 
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
                                CNCu(n1a,n2a,n3a,n4a)=TC(n1b,n1a)*TC(n2b,n2a)*TC(n3b,n3a)*TC(n4b,n4a)*CCu(n1b,n2b,n3b,n4b)+CNCu(n1a,n2a,n3a,n4a); 
                            end
                        end
                    end
                end
            end
        end
    end
end

%do rotations on elascity tensors

CNNb=zeros(3,3,3,3); 
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
                                CNNb(n1a,n2a,n3a,n4a)=TN(n1b,n1a)*TN(n2b,n2a)*TN(n3b,n3a)*TN(n4b,n4a)*CNb(n1b,n2b,n3b,n4b)+CNNb(n1a,n2a,n3a,n4a); 
                            end
                        end
                    end
                end
            end
        end
    end
end

%now we want to calculate and display the 6 form of each elascity tensor

%6 by 6 form
sixNCu=zeros(6,6); 
sixCu=zeros(6,6); 

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
                sixCu(ind1,ind2)=CCu(n1,n2,n3,n4); 
                sixNCu(ind1,ind2)=CNCu(n1,n2,n3,n4); 
            end
        end
    end
end

disp('6 by 6 form of the elascity tensor for Cu'); 
disp(sixNCu); 

%6 by 6 form
sixNNb=zeros(6,6); 
sixNb=zeros(6,6); 

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
                sixNb(ind1,ind2)=CNb(n1,n2,n3,n4); 
                sixNNb(ind1,ind2)=CNNb(n1,n2,n3,n4); 
            end
        end
    end
end

disp('6 by 6 form of the elascity tensor for Nb'); 
disp(sixNNb);


%we need to calculate the P matrices used for enforcing the strain and
%force constraints on the system

P1=zeros(6,6);
P1(1,1)=1; 
P1(2,2)=1; %these give strain constraints
P1(6,6)=1; 
%these give force constraints
for n1=3:1:5
    for n2=1:1:6
        P1(n1,n2)=sixNCu(n1,n2); 
    end
end

P2=zeros(6,6); 
P2(1,1)=1; 
P2(2,2)=1; %these give strain constraints
P2(6,6)=1; 
%these give force constraints
for n1=3:1:5
    for n2=1:1:6
        P2(n1,n2)=sixNNb(n1,n2); 
    end
end

ida=zeros(6,6); 
for n=1:1:6
    ida(n,n)=1; 
end

%display matching matrices
disp('P1'); 
disp(P1); 
disp('P2'); 
disp(P2); 

invP1=ida/P1; 
M=invP1*P2; %use matlab shorthand for inverse matrices

f1=0.5; 
f2=0.5; 

sixComposite=(f1*sixNCu*M+f2*sixNNb)/(f1*M+f2*ida); 

disp('Upper bound on system elascity tensor '); 
disp(sixComposite); 







   
