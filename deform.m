
close all; 
clear all; 
clc; 

%trial 1
%theta=0; 
%phi=0; 

%trial 2
%theta=pi/4; 
%phi=0; 

%trial 3
theta=pi/3; 
phi=pi/3; 

T=[ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ;
    -sin(theta),cos(theta)*cos(phi), sin(phi)*cos(theta) ; 
       0         , -sin(phi) , cos(phi)                 ;   ]; 
   
%we want to calculate the strain field requires to create the deformation
%gradient in this rotated system

strain1=zeros(3,3); 
strain1(2,2)=1; 

strainp=zeros(3,3); 
for n1a=1:1:3
    for n2a=1:1:3
        for n1b=1:1:3
            for n2b=1:1:3
                strainp(n1a,n2a)=T(n1a,n1b)*T(n2a,n2b)*strain1(n1b,n2b)+strainp(n1a,n2a);
            end
        end
    end
end

disp('strain field for rotated system'); 
disp(strainp); 

fp=fopen('in.Cu_deform','w');


fprintf(fp,'echo           log \n'); 
fprintf(fp,'dimension	3 \n'); 
fprintf(fp,'boundary        p p p \n');
fprintf(fp,'units           metal \n'); 
fprintf(fp,'atom_style	atomic \n'); 
fprintf(fp,'neighbor        3.0 bin \n');
fprintf(fp,'neigh_modify    every 1 delay 0 check yes \n');

fprintf(fp,'lattice         fcc 3.615019 \n');
fprintf(fp,'region          box prism 0.0 5.0 0.0 5.0 0.0 5.0 0.0 0.0 0.0 units lattice \n');
fprintf(fp,'create_box      1 box \n');
fprintf(fp,'create_atoms    1 box basis 1 1 \n');

fprintf(fp,'# specify the potential \n');
fprintf(fp,'pair_style      eam/alloy \n');
fprintf(fp,'pair_coeff      * * /home/sskirlo/potentials/CuNbHe_ZBL.eam.alloy Cu #atom 1=Cu, atom 2=Nb, order matters! \n');

fprintf(fp,'# Specify the thermo sampling size/style and output format \n');
fprintf(fp,'thermo          1 \n');
% output full stress and strain vectors
fprintf(fp,'thermo_style    custom step pe xhi yhi zhi xy yz xz pxx pyy pzz pxy pyz pxz \n');

%we need to deform box in specified direction and print out energy after
%each deformation step

%do a 1% deformation along some specified strain vector

%convert strain tensor into vector
st=zeros(1,6); 
st(1,1)=strainp(1,1); 
st(1,2)=strainp(2,2); 
st(1,3)=strainp(3,3); 
st(1,4)=strainp(1,2); 
st(1,5)=strainp(2,3); 
st(1,6)=strainp(1,3); 

dim=[1,1,1]*5;  

%normalize to 1% total strain
st=st*0.001; 

%preforms one deformation step per timestep
fprintf(fp,'fix def all deform 1 x final 0.0 %f y final 0.0 %f z final 0.0 %f xy final %f yz final %f xz final %f units lattice \n',(st(1)+1)*dim(1),(st(2)+1)*dim(2),(st(3)+1)*dim(3),2*st(4)*dim(3),2*st(5)*dim(1),2*st(6)*dim(2)); 
fprintf(fp,'run 1000 \n'); 

%deform in opposite direction
st=-st; 

fprintf(fp,'fix def all deform 1 x final 0.0 %f y final 0.0 %f z final 0.0 %f xy final %f yz final %f xz final %f units lattice \n',(st(1)+1)*dim(1),(st(2)+1)*dim(2),(st(3)+1)*dim(3),2*st(4)*dim(3),2*st(5)*dim(1),2*st(6)*dim(2)); 
fprintf(fp,'run 1000 \n'); 

fclose(fp); 
