
close all; 
clear all; 
clc; 

fp=fopen('in_vac.txt','w'); 

%reads in previous restart file, prevKS1, manipulated by bash script
%open restart file
fprintf(fp,'units           metal \n');
fprintf(fp,'read_restart    prevKS1.restart \n');
fprintf(fp,'lattice         fcc 3.615 \n');
fprintf(fp,'pair_style	eam/alloy \n');
fprintf(fp,'pair_coeff	* * /home/sskirlo/potentials/CuNbHe_ZBL.eam.alloy Cu Nb He #atom 1=Cu, atom 2=Nb, order matters! \n');
fprintf(fp,'thermo          1000 \n');
fprintf(fp,'thermo_style    custom step pe ke etotal temp vol press pxx pyy pzz \n');
fprintf(fp,'thermo 1 \n');
fprintf(fp,'compute kea all ke/atom \n');
fprintf(fp,'compute pea all pe/atom \n');
fprintf(fp,'compute stress all stress/atom \n');

%need to remove atom id specify here, make sure use keyword so doesn't
%shift numbering
fprintf(fp,'group vac id == NA \n'); %id will be specified by bash script 
fprintf(fp,'delete_atoms group vac compress no \n'); %compress keyword set to no means no renumbering 


        
%setup for 2nd minimization step
      
%define handles for region before we minimize z components with atoms inbetween
fprintf(fp,'region h1 block -0.32 97.15 0.124 97.279 30  80  units box \n'); 
fprintf(fp,'region h2 block -0.32 97.15 0.124 97.279 -80 -36 units box \n'); 
fprintf(fp,'group handle1 region h1 \n'); 
fprintf(fp,'group handle2 region h2 \n'); 
fprintf(fp,'fix fp1 all aveforce 0.0 0.0 0.0 region h1 \n'); 
fprintf(fp,'fix fp2 all aveforce 0.0 0.0 0.0 region h2 \n'); 
%leave z component force alone and set others to zero, except in handle region
fprintf(fp,'fix fp all setforce 0.0 0.0 NULL \n'); 

fprintf(fp,'minimize 1.0E-09 1.0E-04 2000 20000 \n'); 

fprintf(fp,'\n \n #do remaining minimizations \n'); 
%now do remaining minimizations 

% do constant volume and then variable volume minimizations
fprintf(fp,'unfix fp \n'); %unfix all atoms zero, so can do minimization with just handles constrained        
fprintf(fp,'minimize 1.0E-08 1.0E-04 2000 20000 \n'); 
fprintf(fp,'unfix fp1 \n');        
fprintf(fp,'unfix fp2 \n');    %unfix handles so can do minimization of cell size
fprintf(fp,'fix bob2 all box/relax x 0.0 y 0.0 \n'); 
fprintf(fp,'minimize 1.0E-08 1.0E-04 2000 20000 \n'); 
fprintf(fp,'unfix bob2 \n \n');  %remove zero pressure constraint

%do annealing and cooling
fprintf(fp,'\n \n #do annealing \n'); 

fprintf(fp,'velocity all create 0.0 3052 dist gaussian \n'); 
fprintf(fp,'velocity all zero linear \n');
fprintf(fp,'fix joe all dt/reset 10 0.000001 0.001 0.005 units lattice #can assume 1 fs timestep \n'); 
%can add in "couple xy" command before drag keyword to make x y dimensions
%change in lockstep
fprintf(fp,'fix bob all nph x 0.0 0.0 0.1 y 0.0 0.0 0.1 xy 0.0 0.0 0.1 drag 2.5 \n'); %allow boundaries to relax

T=0; 
Tstep=20;
fprintf(fp,'\n'); 
for n=1:1:(300)/Tstep
     fprintf(fp,'run 500 pre no post no every 25 "velocity all scale %f" \n ',n*Tstep+T);
end
fprintf(fp,'thermo       10000 \n'); %print out 100 lines during anneal, can moniter
fprintf(fp,'run 10000 pre no post no every 50 "velocity all scale 300" #anneal for 1 ps \n '); %anneal material for 1 ps
Tstep=5; 
for n=1:1:(300)/Tstep
     fprintf(fp,'run 500 pre no post no every 25 "velocity all scale %f" \n ',300-n*Tstep); 
end
fprintf(fp,'run 2000 pre no post no every 25 "velocity all scale %f" \n ',0.0); %stop at 0 K, so know what temperature its at
%unfix nph fix
fprintf(fp,'unfix bob \n');
fprintf(fp,'unfix joe \n'); %unfix time fix also, so can reset timestep counter

fprintf(fp,'\n \n #do final minimization \n'); 

%do final fixed volume minimization
fprintf(fp,'fix fp1 all aveforce 0.0 0.0 0.0 region h1 \n'); 
fprintf(fp,'fix fp2 all aveforce 0.0 0.0 0.0 region h2 \n'); 
fprintf(fp,'minimize 1.0E-10 1.0E-04 2000 20000 \n'); 
fprintf(fp,'unfix fp1 \n'); 
fprintf(fp,'unfix fp2 \n '); 


%write restart file and dump file
fprintf(fp,'write_restart nextKS1.restart \n');  %outputs next.restart, depending on energy state, decides whether uses this restart file or prev
%this needs to be specified by the bash script
fprintf(fp,'reset_timestep NA \n'); 
fprintf(fp,'dump ohsnap all cfg 1 cunb_KS1vac*.dump id type xs ys zs id c_kea c_pea c_stress[1] c_stress[2] c_stress[3] c_stress[4] c_stress[5] c_stress[6] \n');
fprintf(fp,'dump_modify ohsnap element Cu Nb He \n');
%too much to output all dump file images, but can leave space where easy to
%insert
fprintf(fp,'# insert run 1 here  \n');
fprintf(fp,'clear \n');

fclose(fp); 

