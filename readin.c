#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//reads in Demkowicz data file and produces data file lammps can use
int main()
{
	int size=24816; 

	FILE * fp; 
	fp=fopen("CuNb_KS2demk.data","a"); 
	if(fp==NULL)
	{
		printf("bad filename \n"); 
		return 0; 
	}
	FILE * fp2; 
	fp2=fopen("coordinates.txt","r"); 
	if(fp2==NULL)
	{
		printf("bad filename \n"); 
		return 0; 
	}

	//need to write id atom type, and then coordinates
	int n,n1; 
	
	int * types; 
	float * pos; 

	pos=(float *)calloc(size*3,sizeof(float));
	types=(int *)calloc(size,sizeof(int)); 

	float xmin[3]; 
	float xmax[3]; 

	float dimlow[3]; 
	float dimhigh[3]; 

	dimlow[0]=0.0;
	dimlow[1]=0.0;
	dimlow[2]=0.0; 

	dimhigh[0]=97.14494839;
	dimhigh[1]=96.92880912;
	dimhigh[2]=38.05181387;

	for(n=0; n<3; n++)
	{
		xmin[n]=100000;
		xmax[n]=-100000; 
	}

	int swit=-1; //determine index when changes atom type so can print out ids in correct order
	for(n=0; n<size; n++)
	{
		fscanf(fp2,"%d %f %f %f \n",&types[n],&pos[n],&pos[size+n],&pos[2*size+n]);
		if(types[n]==1 && swit==-1)
			swit=n; //will switch only once at index we want

		//we need to wrap around coordinates so the file format is compatible with lammps
		for(n1=0; n1<3; n1++)
		{
			if(pos[n+size*n1]<dimlow[n1])
			{
				pos[n+size*n1]=pos[n+size*n1]+dimhigh[n1]; 
			}
			if(pos[n+size*n1]>dimhigh[n1])
			{
				pos[n+size*n1]=pos[n+size*n1]-dimhigh[n1]; 
			}
		}
		
		for(n1=0; n1<3; n1++)
		{
			if(pos[n+size*n1]<xmin[n1])
			{
				xmin[n1]=pos[n+size*n1]; 
			}
			if(pos[n+size*n1]>xmax[n1])
			{
				xmax[n1]=pos[n+size*n1]; 
			}
		}
		 
	}
	//print out atoms with type 1
	int count=1; 
	for(n=swit; n<size; n++)
	{
		fprintf(fp,"%d %d %f %f %f \n",count,types[n],pos[n],pos[size+n],pos[2*size+n]);
		count++; 
	}
	//print out atoms with type 2
	for(n=0; n<swit; n++)
	{
		fprintf(fp,"%d %d %f %f %f \n",count,types[n],pos[n],pos[size+n],pos[2*size+n]);
		count++;  
	}
	
	printf("min and max dimensions \n"); 
	for(n=0; n<3; n++)
	{
		printf("component %d %f %f \n",n,xmin[n],xmax[n]); 
	}	

	fclose(fp); 
	fclose(fp2); 
	return 0; 
}
