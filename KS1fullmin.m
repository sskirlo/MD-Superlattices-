
close all; 
clear all; 
clc; 

%this is for full minimization of KS1

%dimensions for small simulation
str='s2'; 
xlo=0.0;
xhi=97.2839;
ylo=0.0;
yhi=97.1507; 
zlo=0.0;
zhi=31.3675; 
divide=13.3;




%{
%dimensions for medium simulation
str='m'; 
xlo=0.0;
xhi=97.2839;
ylo=0.0;
yhi=97.1507; 
zlo=0.0;
zhi=163.456; 

%}




a=sprintf('in.CuNb_KS1fullmin%s',str); 

fp=fopen(a,'w'); 

%write input script, start by reading in data file and writing restart file

xrange=0.0; 
yrange=0.0; 

fprintf(fp,'echo           log \n'); 
fprintf(fp,'dimension	3 \n'); 
fprintf(fp,'boundary        p p p \n');
fprintf(fp,'units           metal \n'); 
fprintf(fp,'atom_style	atomic \n'); 
fprintf(fp,'neighbor        3.0 bin \n');
fprintf(fp,'neigh_modify    every 1 delay 0 check yes \n');
a=sprintf('read_data       CuNb_KS1p%s.data \n',str);
fprintf(fp,a);     
a=sprintf('write_restart   cunbks1p%s.restart \n',str); 
fprintf(fp,a);
fprintf(fp,'clear  #clear out all previous input \n'); 

%conduct original gamma-surface minimization
fprintf(fp,' \n \n \n'); 
        
%open restart file
fprintf(fp,'units           metal \n');
a=sprintf('read_restart    cunbks1p%s.restart \n',str); 
fprintf(fp,a);
fprintf(fp,'lattice         fcc 3.615 \n');
fprintf(fp,'pair_style	eam/alloy \n');
fprintf(fp,'pair_coeff	* * /home/sskirlo/potentials/CuNbHe_ZBL.eam.alloy Cu Nb He #atom 1=Cu, atom 2=Nb, order matters! \n');
fprintf(fp,'thermo          1000 \n');
fprintf(fp,'thermo_style    custom step pe ke etotal temp vol press pxx pyy pzz \n');
fprintf(fp,'thermo 1 \n');
fprintf(fp,'compute kea all ke/atom \n');
fprintf(fp,'compute pea all pe/atom \n');
fprintf(fp,'compute stress all stress/atom \n');
        

    fprintf(fp,'\n \n #we make 2 rigid blocks and allow these to shift into equlibrium positions \n \n'); 
    %displace to new position on gamma surface
    %define each of the regions as rigid blocks which can move independently 
    fprintf(fp,'region cu block %f %f %f %f %f %f units box \n',xlo,xhi,ylo,yhi,divide,zhi); 
    fprintf(fp,'region nb block %f %f %f %f %f %f units box \n',xlo,xhi,ylo,yhi,0,divide);
    %create groups of atoms based on original regions (define groups from regions b/c atoms can move outside these bounds
    fprintf(fp,'group copper region cu \n'); 
    fprintf(fp,'group nobium region nb \n'); 
    %fix blocks to only allow rigid translation
    fprintf(fp,'fix fp1 copper aveforce 0.0 0.0 0.0 \n'); 
    fprintf(fp,'fix fp2 nobium aveforce 0.0 0.0 0.0 \n'); 
    %minimize with constraints for all blocks
    fprintf(fp,'minimize 0 0 2000 20000 \n');
    %unfix rigid blocks and do remaining minimizations
    fprintf(fp,'unfix fp1 \n'); 
    fprintf(fp,'unfix fp2 \n');

    fprintf(fp,'\n \n # we do constant volume minization \n \n');
    fprintf(fp,'minimize 0 0 2000 20000 \n'); 


    fprintf(fp,'\n \n # we do variable volume minization \n \n');
    fprintf(fp,'fix bob2 all box/relax x 0.0 y 0.0 z 0.0 \n'); 
    fprintf(fp,'minimize 0 0 2000 20000 \n'); 
    fprintf(fp,'unfix bob2 \n \n');  %remove zero pressure constraint

%do annealing and cooling
fprintf(fp,'\n \n #do annealing \n \n'); 
fprintf(fp,'velocity all create 0.0 3052 dist gaussian \n'); 
fprintf(fp,'velocity all zero linear \n');
fprintf(fp,'fix joe all dt/reset 10 0.000001 0.001 0.005 units lattice #can assume 1 fs timestep \n'); 
fprintf(fp,'fix bob all nph x 0.0 0.0 0.1 y 0.0 0.0 0.1 z 0.0 0.0 0.1 drag 2.5 \n'); %allow boundaries to relax
T=0; 
Tstep=20;
fprintf(fp,'\n');

maxtemp=600; 

for n=1:1:(maxtemp)/Tstep
     fprintf(fp,'run 1000 pre no post no every 25 "velocity all scale %f" \n ',n*Tstep+T);
end
fprintf(fp,'thermo       10000 \n'); %print out 100 lines during anneal, can moniter
fprintf(fp,'run 10000 pre no post no every 50 "velocity all scale 300" #anneal for 1 ps \n '); %anneal material for 1 ps
Tstep=5; 
for n=1:1:(maxtemp)/Tstep
     fprintf(fp,'run 1000 pre no post no every 25 "velocity all scale %f" \n ',maxtemp-n*Tstep); 
end
fprintf(fp,'run 2000 pre no post no every 25 "velocity all scale %f" \n ',0.0); %stop at 0 K, so know what temperature its at
%unfix nph fix
fprintf(fp,'unfix bob \n');
fprintf(fp,'unfix joe \n'); %unfix time fix also, so can reset timestep counter



fprintf(fp,'\n \n # we do variable volume minization \n \n');
fprintf(fp,'fix bob2 all box/relax x 0.0 y 0.0 z 0.0 \n'); 
fprintf(fp,'minimize 0 0 2000 20000 \n'); 
fprintf(fp,'unfix bob2 \n \n');  %remove zero pressure constraint



fprintf(fp,'\n \n # we do final fixed volume minimization \n \n');
fprintf(fp,'minimize 0 0 2000 20000 \n'); 
     

%write restart file and dump file

a=sprintf('write_restart KS1pmined%s.restart \n',str);
fprintf(fp,a);  
fprintf(fp,'reset_timestep 1 \n'); 
a=sprintf('dump ohsnap all cfg 1 cunb_KS1pmined%s*.dump id type xs ys zs id c_kea c_pea c_stress[1] c_stress[2] c_stress[3] c_stress[4] c_stress[5] c_stress[6] \n',str);
fprintf(fp,a);  
fprintf(fp,'dump_modify ohsnap element Cu Nb He \n');
fprintf(fp,'run 1  \n');
fprintf(fp,'clear \n');

fclose(fp); 