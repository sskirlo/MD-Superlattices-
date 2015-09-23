

close all; 
clear all; 
clc; 

fp=fopen('in.CuNb_annealer','w');


fprintf(fp,'units           metal \n');
fprintf(fp,'read_restart    KS1pmined.restart \n');
fprintf(fp,'lattice         fcc 3.615 \n');
fprintf(fp,'pair_style	eam/alloy \n');
fprintf(fp,'pair_coeff	* * /home/sskirlo/potentials/CuNbHe_ZBL.eam.alloy Cu Nb He #atom 1=Cu, atom 2=Nb, order matters! \n');
fprintf(fp,'thermo          1000 \n');
fprintf(fp,'thermo_style    custom step pe ke temp lx ly lz xy yz xz pxx pyy pzz pxy pyz pxz \n');
fprintf(fp,'thermo 1 \n');
fprintf(fp,'compute kea all ke/atom \n');
fprintf(fp,'compute pea all pe/atom \n');
fprintf(fp,'compute stress all stress/atom \n');



fprintf(fp,'velocity all create 0.0 3052 dist gaussian \n'); 
fprintf(fp,'velocity all zero linear \n');
%need to allow expansion in x,y,z direction while maintaining shear
%deformation, will allow us to moniter shear stress

fprintf(fp,'fix bob all nph x 0.0 0.0 0.1 y 0.0 0.0 0.1 z 0.0 0.0 0.1 drag 2.5 \n'); %allow boundaries to relax

Tlevel=10; 
Tstep=2; 
T=0.1; 
Totalsteps=10; 

incs=0; 

fprintf(fp,'run 4000 pre no post no every 50 "velocity all scale %f" \n ',T); 
fprintf(fp,'run 1 \n');
incs=incs+1; 
fprintf(fp,'write_restart anneal%d.restart \n',round(T)); 

%{
for n2=1:1:Totalsteps
    fprintf(fp,'\n'); 
  
    for n=1:1:Tlevel/Tstep
        T=T+Tstep; 
        fprintf(fp,'run 2000 pre no post no every 25 "velocity all scale %f" \n ',T);
        fprintf(fp,'run 1 \n'); %these runs are so we can get stress information
        incs=incs+1; 
    end
    %holding steps
    for n=1:1:Tlevel/Tstep
        fprintf(fp,'run 2000 pre no post no every 50 "velocity all scale %f" \n ',T); 
        fprintf(fp,'run 1 \n');
        incs=incs+1; 
    end
    
    fprintf(fp,'write_restart anneal%d.restart \n',round(T)); %save restart file every 10 K for later deformation
    
end
%}

