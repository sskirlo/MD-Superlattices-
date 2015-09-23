
close all; 
clear all; 
clc; 

xlo=0.0;
xhi=96.8822;
ylo=0.0;
yhi=93.8294; 
zlo=0.0;
zhi=101.716; 
middle=50.0; %coordinate for gap between crystals

fp=fopen('in.CuNb_SPDgsurf','w'); 

%we want to search x y grid of unique displacements and map gamma surface
a0=3.615;  %Cu lattice size

res=0.05; 
xrange=-0:res:2; 
yrange=-0:res:2;
xrange=xrange/a0; 
yrange=yrange/a0; 


%write input script, start by reading in data file and writing restart file

fprintf(fp,'echo           log \n'); 
fprintf(fp,'dimension	3 \n'); 
fprintf(fp,'boundary        p p p \n');
fprintf(fp,'units           metal \n'); 
fprintf(fp,'atom_style	atomic \n'); 
fprintf(fp,'neighbor        3.0 bin \n');
fprintf(fp,'neigh_modify    every 1 delay 0 check yes \n');
fprintf(fp,'read_data       CuNb_SPDps.data \n');     
fprintf(fp,'write_restart   cunbspdgamma.restart \n');
fprintf(fp,'clear  #clear out all previous input \n'); 


sx=size(xrange); 
sx=sx(1,2); 
sy=size(yrange); 
sy=sy(1,2); 

step=0; 
for n1=1:1:sx
    for n2=1:1:sy
        
        fprintf(fp,' \n \n \n'); 
        
        %open restart file
        fprintf(fp,'units           metal \n');
        fprintf(fp,'read_restart    cunbspdgamma.restart \n');
        fprintf(fp,'lattice         fcc 3.615 \n');
        fprintf(fp,'pair_style	eam/alloy \n');
        fprintf(fp,'pair_coeff	* * /home/sskirlo/potentials/CuNbHe_ZBL.eam.alloy Cu Nb He #atom 1=Cu, atom 2=Nb, order matters! \n');
        fprintf(fp,'thermo          1000 \n');
        fprintf(fp,'thermo_style    custom step pe ke etotal temp vol press pxx pyy pzz \n');
        fprintf(fp,'thermo 1 \n');
        fprintf(fp,'compute kea all ke/atom \n');
        fprintf(fp,'compute pea all pe/atom \n');
        fprintf(fp,'compute stress all stress/atom \n');
        
        %setup for first minimization step
        
        %displace to new position on gamma surface
        fprintf(fp,'region co block %f %f %f %f %f %f units box \n',xlo,xhi,ylo,yhi,middle,zhi); 
        fprintf(fp,'region nb block %f %f %f %f %f %f units box \n',xlo,xhi,ylo,yhi,0,middle);
        %create groups of atoms based on original regions (define groups from regions b/c atoms can move outside these bounds)
        fprintf(fp,'group copper region co \n'); 
        fprintf(fp,'group nobium region nb \n'); 
        fprintf(fp,'displace_atoms copper move %f %f 0.0 \n',xrange(n1),yrange(n2));
        %fix blocks to only allow rigid translations in z direction
        fprintf(fp,'fix fp1 copper aveforce NULL NULL 0.0 \n'); 
        fprintf(fp,'fix fp2 nobium aveforce NULL NULL 0.0 \n'); 
        %leave z component force alone and set others to zero
        fprintf(fp,'fix fp all setforce 0.0 0.0 NULL \n'); 
        %now minimize only allowing upper and lower crystal to translate in
        %z direction
        fprintf(fp,'minimize 1.0E-07 1.0E-04 2000 20000 \n'); 
        fprintf(fp,'unfix fp1 \n'); 
        fprintf(fp,'unfix fp2 \n'); 
        fprintf(fp,'unfix fp \n'); 
        
        %fprintf(fp,'minimize 0 0 2000 20000 \n'); 
        
        fprintf(fp,'thermo 1 \n');
        fprintf(fp,'# hello %d \n',step);  %use this as marker 
        fprintf(fp,'run 1  \n');
        fprintf(fp,'clear \n');
        step=step+1; 
    end
end

disp('steps'); 
disp(step); 

fclose(fp); 