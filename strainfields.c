
#define X 0
#define Y 1
#define Z 2
#define ID 3
#define KE 4
#define PE 5
#define SX 6
#define SY 7
#define SZ 8
#define TXY 9
#define TYZ 10
#define TXZ 11
#define TYPE 12  //type field is not output 

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

struct entry
{
	double ent[13]; 
};  //dont forget this 
int printentry(FILE * fp,struct entry * bob,double * dim, double * shear) //will either write to file or print out to screen
{
	char * form; 
	form="%12.6f %12.6f %12.6f %d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n";
	if(fp==NULL)
	{
 		printf(form,bob[0].ent[0],bob[0].ent[1],bob[0].ent[2],(int)bob[0].ent[3],bob[0].ent[4],bob[0].ent[5],bob[0].ent[6],bob[0].ent[7],bob[0].ent[8],bob[0].ent[9],bob[0].ent[10],bob[0].ent[11]);
	}
	//we need to convert coordinates into 0 to 1 form for outputting in dump file
	else
	{
		double x; 
		double y; 
		double z;
		z=bob[0].ent[2]/dim[2]; 
		y=bob[0].ent[1]-bob[0].ent[2]*shear[1]/dim[2]; 
		x=bob[0].ent[0]-bob[0].ent[2]*shear[2]/dim[2]-y*shear[0]/dim[1];
		x=x/dim[0];
		y=y/dim[1]; 
		if(x>1.0 || x<0.0 || y>1.0 || y<0.0 || z>1.0 || z<0.0 )
		{
  			printf("bad atom coordinate exceeds range id %d \n",(int)bob[0].ent[ID]); 
		}
    	fprintf(fp,form,x,y,z,(int)bob[0].ent[3],bob[0].ent[4],bob[0].ent[5],bob[0].ent[6],bob[0].ent[7],bob[0].ent[8],bob[0].ent[9],bob[0].ent[10],bob[0].ent[11]);
	}
	return 0;  
}
int compareid(const void *ptr1, const void *ptr2) //compare function to be passed to qsort
{
	double e1, e2, an1, an2; 
	struct entry* p1;
	struct entry* p2; 
	p1=(struct entry*)ptr1; 
	p2=(struct entry*)ptr2;  
	e1=p1[0].ent[ID]; 
	e2=p2[0].ent[ID];
	if(e1>e2)
		return 1;
	else if(e1<e2)
		return -1;
	else
		return 0;
}
int comparetype(const void *ptr1, const void *ptr2) //compare function to be passed to qsort
{
	double e1, e2, an1, an2; 
	struct entry* p1;
	struct entry* p2; 
	p1=(struct entry*)ptr1; 
	p2=(struct entry*)ptr2;  
	e1=p1[0].ent[TYPE]; 
	e2=p2[0].ent[TYPE];
	if(e1>e2)
		return 1;
	else if(e1<e2)
		return -1;
	else
		return 0;
}
//code works by checking that ids are properly ordered and that there are no duplicates, some ids are deleted
int checkAtoms(struct entry * dat, int atomnumber)// purpose of command is to check that all atom ids present
{
	int count;
	int prev; 
	prev=0;
	for(count=0; count<atomnumber; count++)
	{
		//printf("%d \n",(int)dat[count].ent[ID]); 
		if( ((int)(dat[count].ent[ID]))<=prev ) //if not ascending at least one, something went wrong
		{
			printentry(NULL,&dat[count],NULL,NULL); 
			printf("missing atom %d \n", (count+1)); 
			return 1;
		} 
		prev=dat[count].ent[ID];
	}

	printf("have %d of %d atoms accounted for \n", count,atomnumber); 
	return 0; 	
}
//calculate distance don't take into account boundaries
double distance(struct entry b1, struct entry b2)
{
	double d; 
	d=pow(b1.ent[0]-b2.ent[0],2.0)+pow(b1.ent[1]-b2.ent[1],2.0)+pow(b1.ent[2]-b2.ent[2],2.0);
	d=pow(d,0.5); 
	return d; 
}
//read in data from file fp and label atom type 1st velocity field
int readin(FILE *fp,int size, struct entry * buf, double * dim, double * shear)
{
	char * formin; 
	formin="%f %f %f %f %f %f %f %f %f %f %f %f \n"; 
	
	int i,n,check; 
	float im[12]; 
	
	int atomtype=1;
	int changed=0;  
	
	for(n=0;n<size;n++)
	{
		check=fscanf(fp,formin,&im[0],&im[1],&im[2],&im[3],&im[4],&im[5],&im[6],&im[7],&im[8],&im[9],&im[10],&im[11]); 
		if(check!=12)
		{
			changed=1; 
			fgets(NULL,0,fp); 
			fgets(NULL,0,fp); //go to next line	
			fscanf(fp,formin,&im[0],&im[1],&im[2],&im[3],&im[4],&im[5],&im[6],&im[7],&im[8],&im[9],&im[10],&im[11]);
			atomtype=2; //change atom type
			printf("changed atom type at %d of %d atoms \n",n,size); 
			printf("id of changed atom %f \n",im[ID]); 
		}
		for(i=0; i<12; i++)
		{
			buf[n].ent[i]=(double)im[i]; //absolute addressing based on ID
		}
		for(i=0; i<3; i++) //positions always here
		{
			buf[n].ent[i]=buf[n].ent[i]*dim[i]; //convert from relative to absolute coordinates
		}
		//we need to add shear components to each atom from reading in box coordinates
		buf[n].ent[0]=buf[n].ent[1]/dim[1]*shear[0]+buf[n].ent[2]/dim[2]*shear[2]+buf[n].ent[0]; 
		buf[n].ent[1]=buf[n].ent[2]/dim[2]*shear[1]+buf[n].ent[1]; 
   
		buf[n].ent[TYPE]=atomtype; //put in atom type

		/*if(changed==1)
			printentry(NULL,&buf[n]);
		changed=0;*/ 
	}

	return 0; 
}
//this method writes out full dumpfile
int writeoutall(char * str, struct entry * d, int count, double * dim, double * shear, char * app)
{
	int n; 
	char text[100]; 
	int size=100; 
	//need to first open up original input file so can copy over key parts of formatting
	FILE * fp; 
	fp=fopen(str,"r");
	char mat[50]; 
	join(mat,"dumpy",app,".dump");
	FILE * out; 
	out=fopen(mat,"w"); 
	
	fprintf(out,"Number of particles = %d \n",count);
	fgets(text,100,fp);  		 
	for(n=0; n<23; n++)
	{
		fgets(text,100,fp); 
		fprintf(out,"%s",text); 	
	}
	printf("sorting \n"); 		
	qsort((void*)d,count,sizeof(struct entry),comparetype); //need to compare atom ids by type
	writeout(out,d,count,dim,shear); 
}
//method will write out 1st atom type, and if another atom type appears will automatically begin writing second part
//write-out, write-out, write for ruin, and the world's ending
int writeout(FILE *fp, struct entry * d,  int count, double * dim, double * shear)
{
	int n; 
	int flag=0;
	printf("writing out ids \n"); 
	for(n=0;n<count;n++)
	{
		if( (((int)(d[n].ent[TYPE]))==2)&& (flag==0) )//if at second atom type print out appropriate lines
		{
			fprintf(fp, "92.9064 \n"); 
		    	fprintf(fp,"Nb \n"); //get two middle lines one has mass other has atom type 
		    	flag=1; //do this so won't keep printing out
		}
		printentry(fp,&d[n],dim,shear); 
	}
	return 0; 
}

//puts together three strings in master char 
int join(char* mas,char * one, char * two, char *three)
{
	//printf("joining \n");
	//mas is 30 chars long, so will be able to hold all the data
	char* joe[3]; 
	joe[0]=one; 
	joe[1]=two; 
	joe[2]=three;
	int n1,n2;
	int len;
	int pos;  
	pos=0; 
	for(n1=0; n1<3; n1++)
	{
		for(n2=0; n2<strlen(joe[n1]); n2++)
		{
				mas[pos]=*(joe[n1]+n2); 
				pos++; 
		}
	}
	mas[pos]='\0'; 
	//printf("finish joining \n");
	printf("%s \n",mas);  
}

//if this code returns -1, need to terminate program operation b/c compromised
int readinall(char * str ,struct entry ** dat,int * atoms,double * dim, double * shear )
{

	int n1,n; 	
	FILE * fp=fopen(str,"r"); 
	if(fp==NULL)
	{
		printf("Invalid input"); 
		return -1;
	}
	 
	int buff1=100; 
	char text[100]; //generally safe assumption that all lines manipulating are shorter than 100 characters      	
	
	for(n1=0; n1<4; n1++)
		fscanf(fp,"%s",text); 
	
	fscanf(fp," %d \n ",atoms); //-1  //get number of atoms

	printf("%d atoms in system with 2 types \n",(*atoms));
	float tdim;
  
	//start getting dimensions of box
	fgets(text,buff1,fp); //-1 line
	fscanf(fp, "%s",text); 
	fscanf(fp, "%s",text); //get two strings out of line
	fscanf(fp, "%f \n",&tdim); //read line data into dimension -1
	dim[0]=(double)tdim; //need to be careful with readins, can't read in doubles, only floats, so have to convert
	for(n=0; n<3; n++)  
	{
		fgets(text,buff1,fp); //-3 lines
	}
	fscanf(fp, "%s",text); 
	fscanf(fp, "%s",text); //get two strings out of line
	fscanf(fp, "%f \n",&tdim); //read line data into dimension -1
	shear[0]=(double)tdim; //need to be careful with readins, can't read in doubles, only floats, so have to convert
	fgets(text,buff1,fp); 
	fscanf(fp, "%s",text); 
	fscanf(fp, "%s",text); //get two strings out of line
	fscanf(fp, "%f \n",&tdim); //read line data into dimension -1
	dim[1]=(double)tdim; //need to be careful with readins, can't read in doubles, only floats, so have to convert
	fgets(text,buff1,fp);
	fgets(text,buff1,fp);
	fscanf(fp, "%s",text); 
	fscanf(fp, "%s",text); //get two strings out of line
	fscanf(fp, "%f \n",&tdim); //read line data into dimension -1
	shear[2]=(double)tdim; //need to be careful with readins, can't read in doubles, only floats, so have to convert
	fgets(text,buff1,fp);
	fscanf(fp, "%s",text); 
	fscanf(fp, "%s",text); //get two strings out of line
	fscanf(fp, "%f \n",&tdim); //read line data into dimension -1
	shear[1]=(double)tdim; //need to be careful with readins, can't read in doubles, only floats, so have to convert
	fgets(text,buff1,fp); 
	fscanf(fp, "%s",text); 
	fscanf(fp, "%s",text); //get two strings out of line
	fscanf(fp, "%f \n",&tdim); //read line data into dimension -1
	dim[2]=(double)tdim; //need to be careful with readins, can't read in doubles, only floats, so have to convert
	for(n=0; n<4; n++)  
	{
		fgets(text,buff1,fp); //-3 lines
	}
	//end getting dimensions of box

	printf("The dimensions are %f A, %f A, %f A %f A %f A %f A \n", dim[0],dim[1],dim[2],shear[0],shear[1],shear[2]);
	
	for(n=0; n<10; n++)//take out remaining lines
	{
		fgets(text,buff1,fp); 
	}
	
	(*dat)=(struct entry *)calloc(*atoms,sizeof(struct entry)); //get data to store all the atom information
	if((*dat)==NULL)
	{
		printf("not enough data allocated \n"); 
		return -1; 
	}
	
	readin(fp,(*atoms),*dat,dim,shear); //reads in data starting at rd1 
	
	printf("sorting atoms by id \n"); 
	qsort((void*)(*dat),(*atoms),sizeof(struct entry),compareid); //need to sort data by id, so can check that everything read in properly

	printf("data read in \n");  
	n=checkAtoms((*dat),(*atoms)); 
	
	return n; 
}
int minimum(struct entry * d1, int size, int component,double * min)
{
	int n; 
	(*min)=1000000; 
	for(n=0; n<size; n++)
	{
		if( (*min)>d1[n].ent[component])
		{
			(*min)=d1[n].ent[component]; 
		}		
	}
	return 0; 
}
int maximum(struct entry * d1, int size, int component,double * max)
{
	int n; 
	(*max)=-1000000; 
	for(n=0; n<size; n++)
	{
		if( (*max)<d1[n].ent[component])
		{
			(*max)=d1[n].ent[component]; 
		}		
	}
	return 0; 
}
//calculate vector between two points, also returns distance
double vec(struct entry * one, struct entry * two, double * vecs, double * dimr)
{	
	double d[3];
	int n; 
	for(n=0; n<3; n++)
		d[n]=((*one).ent[n]-(*two).ent[n]);
 
	for(n=0; n<3; n++)//if any coordinates more than half the distance apart, know not neighboring cells
	{
		if(d[n]>dimr[n]/2 || d[n]<-dimr[n]/2) //if displacement is more than half the simulation cell length, we know we have loop around
		{
			if(d[n]>0)
			{
				d[n]=d[n]-dimr[n]; 
			}
			else
			{
				d[n]=d[n]+dimr[n]; 
			}
		}	
	}
	for(n=0;n<3;n++)
		vecs[n]=d[n];
	for(n=0;n<3;n++)
		d[n]=pow(d[n],2.0);
	d[0]=pow(d[0]+d[1]+d[2],0.5);
	return d[0];  //return distance
}
//calculate the disregistry 
int calculateDisregistry(struct entry *rd1, int size, double *dmir, double * shear,int type, int includestressfields, char * append)
{
	int n,n2; 
	FILE * fp; 
	FILE * out; 
	char filename[50]; 
	sprintf(filename,"refdisreg%d.txt",type); 
	fp=fopen(filename,"r"); 
	join(filename,"disregistry",append,".txt");
	out=fopen(filename,"w");
	if(fp==NULL || out==NULL)
	{
		printf("not able to open disregistry reference file \n"); 
		return 0; 
	}	 
	int number; 
	fscanf(fp,"%d \n",&number);
	printf("reading in disregistry vectors, number reading in %d \n",number); 
	printf("include stress fields %d \n",includestressfields); 
	
	int id1,id2; 
	double vecsor[3];
	double vecsd0[3]; 	 
	double vecsd[3]; 

	for(n=0; n<number; n++)
	{
		fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf \n",&id1,&id2,&(vecsd0[0]),&(vecsd0[1]),&(vecsd0[2]),&(vecsor[0]),&(vecsor[1]),&(vecsor[2])); 
		//now we need to look up ids and calculate difference and subtract this vector, list is already sorted by id number		
		for(n2=0; n2<3; n2++)
		{
			vecsd[n2]=rd1[id1-1].ent[n2]-rd1[id2-1].ent[n2]-vecsd0[n2];  
		}
		if( (abs(vecsd[0])>dmir[0]/2) || (abs(vecsd[1])>dmir[1]/2) || (abs(vecsd[2])>dmir[2]/2) )
		{
			continue; //we want to not plot a result if the disregistry is more than half the box size
		} 
		else
		{
			if(includestressfields==0)
			{
				fprintf(out,"%d %d %f %f %f %f %f %f \n",(int)rd1[id1-1].ent[ID],(int)rd1[id2-1].ent[ID],vecsd[0],vecsd[1],vecsd[2],vecsor[0],vecsor[1],vecsor[2]);
			}
			if(includestressfields==1)
			{
				fprintf(out,"%d %d %f %f %f %f %f %f %f %f %f %f %f %f \n",(int)rd1[id1-1].ent[ID],(int)rd1[id2-1].ent[ID],vecsd[0],vecsd[1],vecsd[2],vecsor[0],vecsor[1],vecsor[2],rd1[id1-1].ent[SX],rd1[id1-1].ent[SY],rd1[id1-1].ent[SZ],rd1[id1-1].ent[TXY],rd1[id1-1].ent[TYZ],rd1[id1-1].ent[TXZ]); 
			}
		}
	} 	
	
	fclose(out); 
	fclose(fp); 
	return 0; 
}
//need to create a reference file for the disregistry calculation
//this file determines the location of the interface, preforms a deformation on the Cu atoms, and creates a based on nearest neighbor atoms
//the output file will contain the pairs of Cu and Nb ids and half the value of their displacement vectors 

//zhi, zlo specify the width of the regions we are getting the disregistry vectors for
//type specifies whether deforming Cu side of interface (1) or Nb side (2) in calculating disregistry
int createRefFile(struct entry *rd1, int size, double *dmir, double * shear, double zlo, double zhi,int type,int steps, double * steplocations, double stepheight, char * str, char * append)
{
	int n,n1,n2,n3,countnb,countcu,nstep; 
	double acu=3.615019; 
	double anb=3.3008; 
	double add; 
	countnb=0; 
	countcu=0; 
	
	if(type!=4)
	{
		//we want to count the number of each type of atom in the interface
		for(n=0; n<size; n++)
		{
			if(rd1[n].ent[Z]>zlo && rd1[n].ent[Z]<zhi)
			{
				if(rd1[n].ent[TYPE]==1)
					countcu++; 
				else
					countnb++; 	
			}	
		}
	}
	else
	{
		printf("we are reading in interface atoms with step modulation \n"); 
		printf("step number is %d step height is %f \n",steps,stepheight); 
		printf("The following are the locations of the steps: \n"); 
		for(n=0; n<(steps*2); n++)
		{
			printf("%f \n",steplocations[n]); 
		}
	
	//we need to shift zlo and zhi up and down by step height to deal with modulation
		for(n=0; n<size; n++)
		{
			add=0; //as default add is zero, if we find lies in range with shifted terrace, we change add to step size
			for(nstep=0; nstep<steps; nstep++)
			{
				if( (steplocations[nstep*2]<rd1[n].ent[Y]) && (steplocations[nstep*2+1]>rd1[n].ent[Y]))
				{
					add=stepheight; 
				}
			}
			if(rd1[n].ent[Z]>(zlo-add) && rd1[n].ent[Z]<(zhi-add))
			{
				if(rd1[n].ent[TYPE]==1)
					countcu++; 
				else
					countnb++; 	
			}	
		}
	}
	
	
	printf("%d atoms on copper side of interface, %d on Nb side \n",countcu,countnb); 
	
	struct entry * cuatoms; 
	struct entry * nbatoms; 
	cuatoms=(struct entry *)calloc(countcu,sizeof(struct entry)); 
	nbatoms=(struct entry *)calloc(countnb,sizeof(struct entry)); 
	
	//we need to reset counters to 0
	countcu=0; 
	countnb=0; 
	//now we need to fill up the lists with each type of atom
	if(type!=4)
	{
		for(n=0; n<size; n++)
		{
			if(rd1[n].ent[Z]>zlo && rd1[n].ent[Z]<zhi)
			{
				if(rd1[n].ent[TYPE]==1)
				{
					cuatoms[countcu]=rd1[n]; 
					countcu++;
				} 
				else
				{
					nbatoms[countnb]=rd1[n]; 
					countnb++;
				} 	
			}	
		}
	}
	else
	{
		for(n=0; n<size; n++)
		{
			add=0; //as default add is zero, if we find lies in range with shifted terrace, we change add to step size
			for(nstep=0; nstep<steps; nstep++)
			{
				if( (steplocations[nstep*2]<rd1[n].ent[Y]) && (steplocations[nstep*2+1]>rd1[n].ent[Y]))
				{
					add=stepheight; 
				}
			}
			if(rd1[n].ent[Z]>(zlo-add) && rd1[n].ent[Z]<(zhi-add))
			{
				if(rd1[n].ent[TYPE]==1)
				{
					cuatoms[countcu]=rd1[n]; 
					countcu++;
				} 
				else
				{
					nbatoms[countnb]=rd1[n]; 
					countnb++;
				} 	
			}	
		}
	}
	//now we need to preform a transformation on either the Cu side or Nb side to create a coherent interface
	
	double F1[2][2]; 
	F1[0][0]=-(pow(6.0,0.5)-3)*(acu+anb)/acu; 
	F1[0][1]=(pow(3.0,0.5)-pow(2.0,0.5))*(acu-anb)/acu; 
	F1[1][0]=(6-2*pow(6.0,0.5))/pow(2.0,0.5)+3*(2-pow(6.0,0.5))/pow(2.0,0.5)*anb/acu; 
	F1[1][1]=(pow(6.0,0.5)-2)-(pow(6.0,0.5)-3)*anb/acu; 
	double F2[2][2]; 
	F2[0][0]=(pow(6.0,0.5)-2)*(acu+anb)/anb; 
	F2[0][1]=(pow(3.0,0.5)-pow(2.0,0.5))*(acu-2*anb)/anb/2; 
	F2[1][0]=(6-3*pow(6.0,0.5))/pow(3.0,0.5)+2*(3-pow(6.0,0.5))/pow(3.0,0.5)*acu/anb; 
	F2[1][1]=(3-pow(6.0,0.5))+(pow(6.0,0.5)/2-1)*acu/anb; 

        /*
	printf("we calculated transformation compare each element \n"); 
	printf(" %f  \n ",F1[0][0]); 
	printf(" %f  \n ",F1[1][0]); 
	printf(" %f  \n ",F1[0][1]); 
	printf(" %f  \n ",F1[1][1]); 
	printf(" %f  \n ",F2[0][0]); 
	printf(" %f  \n ",F2[1][0]); 
	printf(" %f  \n ",F2[0][1]); 
	printf(" %f  \n ",F2[1][1]); 
	*/

	//now we need to apply the transformation to either the Cu layer or the Nb layer
	 
	double pos[2]; 
	if(type==1)
	{
		printf("type 1 disregistry calculation, applying transformation to Cu side \n"); 
		for(n=0; n<countcu; n++)
		{
			pos[0]=0; 
			pos[1]=0; 
			pos[0]=F1[0][0]*cuatoms[n].ent[X]+pos[0]; 
			pos[0]=F1[0][1]*cuatoms[n].ent[Y]+pos[0]; 
			pos[1]=F1[1][0]*cuatoms[n].ent[X]+pos[1]; 
			pos[1]=F1[1][1]*cuatoms[n].ent[Y]+pos[1]; 

			//printf("initial position %f %f %d id \n",cuatoms[n].ent[X],cuatoms[n].ent[Y],(int)cuatoms[n].ent[ID]); 
			//printf("final position %f %f \n",pos[0],pos[1]); 
			cuatoms[n].ent[X]=pos[0]; 
			cuatoms[n].ent[Y]=pos[1]; 
		}
	}
	//other transformation
	if(type==2)
	{
		printf("type 2 disregistry calculation, applying transformation to Nb side \n"); 
	}
	//no transformation
	if(type==3 || type==4)
	{
		printf("type 3 or type 4 disregistry calculation, no transformation being applied \n"); 	
	}
	
	//now we need to match each nb atom to a cu atom, we will always do nb to cu atoms, because there are more cu atoms than nb atoms
	//the net effect of the transformation will always be to expand cu atom interface relative to nb interface
	
	
	char filename[50]; 
	sprintf(filename,"refdisreg%d.txt",type); 
	printf("writing out disregistry reference %s \n",filename); 

	FILE * fp; 
	fp=fopen(filename,"w"); 
	
	int sid; 
	double dist,min;
	double diff[3]; 
	double origin[3];  
	fprintf(fp,"%d \n",countnb);  
	for(n=0; n<countnb; n++)
	{
		min=1000000; 
		sid=-1; 
		//find minimum distance to current nb atom
		for(n1=0; n1<countcu; n1++)
		{
			dist=distance(nbatoms[n],cuatoms[n1]); 
			if(dist<min)
			{
				min=dist; 
				sid=n1; 
			}	
		}	
		//now we want to print out id pairs and vector between them
	
		for(n1=0; n1<3; n1++)
		{
			diff[n1]=nbatoms[n].ent[n1]-cuatoms[sid].ent[n1]; 
			origin[n1]=(nbatoms[n].ent[n1]+cuatoms[sid].ent[n1])/2; 
		}
		fprintf(fp,"%d %d %f %f %f %f %f %f \n",(int)nbatoms[n].ent[ID],(int)cuatoms[sid].ent[ID],diff[0],diff[1],diff[2],origin[0],origin[1],origin[2]); 
	}
	
	fclose(fp); 

	//we want to output dump file for visualization, we need to combine atoms together into single struct
	//we also need to rescale coordinates so everything will fall within boundaries

	printf("writing out dump file \n"); 
	
	struct entry * outstruct; 
	outstruct=(struct entry *)calloc(countnb+countcu,sizeof(struct entry));
	if(outstruct==NULL)
	{
		printf("invalid memory address \n"); 
		return 0; 
	}
	int count=0;  
	for(n=0; n<countcu; n++)
	{
		if(cuatoms[n].ent[X]>0 && cuatoms[n].ent[X]<dmir[0] && cuatoms[n].ent[Y]>0 && cuatoms[n].ent[Y]<dmir[1])
		{
			outstruct[count]=cuatoms[n]; 
			count++; 
		}
	}
	printf("Cu atoms in boundary %d \n",count); 
	for(n=0; n<countnb; n++)
	{
		if(nbatoms[n].ent[X]>0 && nbatoms[n].ent[X]<dmir[0] && nbatoms[n].ent[Y]>0 && nbatoms[n].ent[Y]<dmir[1])
		{
			outstruct[count]=nbatoms[n]; 
			count++; 
		}
	}
	printf("Total atoms in boundary %d \n",count); 
	writeoutall(str,outstruct,count,dmir,shear,append); 

	return 0; 	
}
//want to calculate the average spacing between layers
int calculateAvgZ(struct entry* d1, int size, double * dmir)
{
	//need to calculate min and max z coordinates before anything else
	int n,n2;
	double zhi=-100000; 
	double zlo=100000; 

	for(n=0; n<size; n++)
	{
		if(zhi<d1[n].ent[Z])
		{
			zhi=d1[n].ent[Z]; 
		}
		if(zlo>d1[n].ent[Z])
		{
			zlo=d1[n].ent[Z]; 
		}
	}
	printf("zhi is %lf and zlo is %lf \n",zhi,zlo); 

	double a01=3.615; 
	double a02=3.3008; 
	double estlayer=(a01/pow(3,0.5)+a02/pow(2,0.5))/2.0; 
	double bindist=1.0; //set binning distance
	int cells=(int)((zhi-zlo)/bindist)+2; 
	
	int * counter; 
	double * accum; 
	
	//need to create accumulator and counter for z distances
	counter=(int *)calloc(cells,sizeof(int)); 
	accum=(double *)calloc(cells,sizeof(double)); 
	int pos; 
	
	//initialize data
	for(n=0; n<cells; n++)
	{
		counter[n]=0; 
		accum[n]=0; 	
	}
	for(n=0; n<size; n++)
	{
		pos=(int)((d1[n].ent[Z]-zlo)/bindist); 
		counter[pos]++; 
		accum[pos]=d1[n].ent[Z]+accum[pos]; 
	}
	
	//we want to display all relevant statistics as we conduct different operations
	/*for(n=0; n<cells; n++)
	{
		double avg=0; 
		if(counter[n]!=0)
		{
			avg=accum[n]/counter[n]; 
		}	
		printf("%d %lf \n",counter[n],avg); 
	}*/
	//we want to remove 0s between groups, first count up number of nonzero counts
	int countzero=0; 
	for(n=0; n<cells; n++)
	{
		if(counter[n]!=0)
		{
			countzero++; 	
		}
	}
	printf("countzero %d \n",countzero); 
	//need to identify current nonzero one and move to start of list, don't worry about overwriting
	for(n=0; n<countzero; n++)
	{
		for(n2=n; n2<cells;n2++)
		{
			if(counter[n2]!=0)
				break; 	
		}
		//printf("n2 %d \n",n2); 
		counter[n]=counter[n2]; //transfer data
		accum[n]=accum[n2];
		if(n!=n2)//otherwise we don't want to set to 0
		{ 
			counter[n2]=0; //erase in old place
			accum[n2]=0;
		} 
	}
	printf("condensed z position data %d countzero \n",countzero); 
	for(n=0; n<countzero; n++)
	{
		double avg=0; 
		if(counter[n]!=0)
		{
			avg=accum[n]/counter[n]; 
		}	
		printf("%d %lf \n",counter[n],avg); 
	}
	//now we need to check if the average between any of the groups is too close, so we can merge results together
	
	/*double tol=0.5; 
	for(n=0; n<countzero-1; n++)
	{
		if((accum[n]/counter[n])>(-tol+accum[n+1]/counter[n+1]))
		{
			accum[n]=accum[n]+accum[n+1]; 
			counter[n]=counter[n+1]+counter[n]; 
			//need to shift everything else down by 1
			for(n2=n+2;n2<countzero;n2++)
			{
				counter[n2-1]=counter[n2]; 
				accum[n2-1]=accum[n2]; 	
			}	
			countzero=countzero-1; //need to update b/c list got shorter b/c of operation
		}
	}*/
	printf("merging cells by minimum size cutoff method \n"); 
	int minimum=1000; 
	//another method is to merge cells automatically if they have fewer than 1000 atoms
	for(n=0; n<=countzero-1; n++)
	{
		if(counter[n]<minimum) //merge with left or right, depending on which average is closer
		{
			double avg=accum[n]/counter[n]; 
			double avgb;
			double avga; 
			if(n==0)  //particular cases to decide which average to pick
			{
				avgb=(accum[countzero-1]/counter[countzero-1])-dmir[2]; 
			}
			else
			{
				avgb=(accum[n-1]/counter[n-1]); 
			}
			if(n==countzero-1)  //particular cases to decide which average to pick
			{
				avga=(accum[0]/counter[0])+dmir[2]; 
			}
			else
			{
				avga=(accum[n+1]/counter[n+1]); 
			}
			//now we need to determine which average is closer
			if((avg-avgb)>(avga-avg))  //merge to the one in front
			{
				if(n==countzero-1)
				{
					accum[0]=accum[0]+accum[n]-dmir[2]*counter[n];  //need to wrap around to other dimension 
					counter[0]=counter[0]+counter[n];
				} 
				else
				{
					accum[n+1]=accum[n+1]+accum[n];   
					counter[n+1]=counter[n+1]+counter[n];
				}
			}
			else  //merge to the one in back
			{
				if(n==0)
				{
					accum[countzero-1]=accum[countzero-1]+accum[n]+dmir[2]*counter[n];  //need to wrap around to other dimension 
					counter[countzero-1]=counter[countzero-1]+counter[n];
				} 
				else
				{
					accum[n-1]=accum[n-1]+accum[n];   
					counter[n-1]=counter[n-1]+counter[n];
				}
			}
			//need to shift everything else down by 1
			for(n2=n+1;n2<countzero;n2++)
			{
				counter[n2-1]=counter[n2]; 
				accum[n2-1]=accum[n2]; 	
			}	
			countzero=countzero-1; //need to update b/c list got shorter b/c of operation
			n=n-1; //need to go back one, b/c shifted new one to the old n position, need to check it too		
		}
	}
	//we can display final version of list
	printf("final rows prepared %d countzero \n",countzero); 
	for(n=0; n<countzero; n++)
	{
		double avg=0; 
		if(counter[n]!=0)
		{
			avg=accum[n]/counter[n]; 
		}	
		printf("%d %lf \n",counter[n],avg); 
	}
	
	double accumavg=0; 
	//we want to calculate the average lattice spacing
	for(n=0; n<countzero-1; n++)
	{
		accumavg=accumavg+(accum[n]/counter[n]-accum[n+1]/counter[n+1]); 
	}
	//last row is special case
	accumavg=accumavg+(accum[countzero-1]/counter[countzero-1]-accum[0]/counter[0]-dmir[2]);  //need to subtract away box dimension 	

	accumavg=accumavg/(countzero); 
	
	printf("calculated average lattice spacing as %lf \n",accumavg); 
	printf("estimated average is %lf \n",estlayer); 
	
	//we also want to calculate the average of the interface spacings and the regular lattice independently
	int thresh=200; //if accum changes by this much, consider interface average

	double accumavgbulk=0; 
	double accumavgint=0;
	int countint=0; 
	int countbulk=0;  
	for(n=0; n<countzero; n++)
	{
		if(n!=(countzero-1))
		{
			if( ((counter[n]-counter[n+1])>thresh) || ((counter[n+1]-counter[n])>thresh) )
			{
				printf("%d \n",counter[n]-counter[n+1]); 
				accumavgint=(accum[n]/counter[n]-accum[n+1]/counter[n+1])+accumavgint;
				countint++; 
			}
			else
			{
				accumavgbulk=(accum[n]/counter[n]-accum[n+1]/counter[n+1])+accumavgbulk;
				countbulk++; 
			}
		}
		else //have special case for endpoints
		{
			if( ((counter[n]-counter[0])>thresh) || ((counter[0]-counter[n])>thresh) )
			{
				accumavgint=(accum[n]/counter[n]-accum[0]/counter[0]-dmir[2])+accumavgint;
				countint++; 
			}
			else
			{
				accumavgbulk=(accum[n]/counter[n]-accum[0]/counter[0]-dmir[2])+accumavgbulk;
				countbulk++; 
			}
		}
	}	

	printf("calculated bulk and interface row spacing averages seperately \n"); 
	printf("number of interfaces %d and number of rows %d \n",countint,countbulk); 
	printf("interface average %lf bulk average %lf \n",accumavgint/countint,accumavgbulk/countbulk); 
	

	//we also want to calculate the average interface spacing by atom type with interface effects removed
	thresh=200; //if accum changes by this much, consider interface average
	int count1=1200; 
	int count2=1600; //consider average lattice constant for individual rows 

	double accumavgbulkcu=0;
	double accumavgbulknb=0;  
	accumavgint=0;
	countint=0; 
	int countbulkcu=0;  
	int countbulknb=0;  
	for(n=0; n<countzero; n++)
	{
		if(n!=(countzero-1))
		{
			if( ((counter[n]-counter[n+1])>thresh) || ((counter[n+1]-counter[n])>thresh) )
			{
				printf("%d \n",counter[n]-counter[n+1]); 
				accumavgint=(accum[n]/counter[n]-accum[n+1]/counter[n+1])+accumavgint;
				countint++; 
			}
			else
			{
				if( ((counter[n]+thresh)>count1) && ((counter[n]-thresh)<count1) )
				{
					accumavgbulkcu=(accum[n]/counter[n]-accum[n+1]/counter[n+1])+accumavgbulkcu;
					countbulkcu++;
				} 
				else
				{
					accumavgbulknb=(accum[n]/counter[n]-accum[n+1]/counter[n+1])+accumavgbulknb;
					countbulknb++;
				}
			}
		}
		else //have special case for endpoints
		{
			if( ((counter[n]-counter[0])>thresh) || ((counter[0]-counter[n])>thresh) )
			{
				accumavgint=(accum[n]/counter[n]-accum[0]/counter[0]-dmir[2])+accumavgint;
				countint++; 
			}
			else
			{
				if( ((counter[n]+thresh)>count1) && ((counter[n]-thresh)<count1) )
				{
					accumavgbulkcu=(accum[n]/counter[n]-accum[0]/counter[0]-dmir[2])+accumavgbulkcu;
					countbulkcu++;
				} 
				else
				{
					accumavgbulknb=(accum[n]/counter[n]-accum[0]/counter[0]-dmir[2])+accumavgbulknb;
					countbulknb++;
				} 
			}
		}
	}	

	printf("calculated bulk cu and nb and interface row spacing averages seperately \n"); 
	printf("number of interfaces %d and number of cu rows %d and nb rows %d \n",countint,countbulkcu,countbulknb); 
	printf("interface average %lf cu bulk average %lf and nb bulk average %lf \n",accumavgint/countint,accumavgbulkcu/countbulkcu,accumavgbulknb/countbulknb); 
	printf("expect nb average %lf and expected cu average %lf \n",a01/pow(3,0.5),a02/pow(2,0.5) ); 

	int sumcheck=0; 
	for(n=0; n<countzero; n++)
	{
		sumcheck=counter[n]+sumcheck; 
	}

	printf("atom number is conserved %d, %d \n",sumcheck,size); 

	return 0; 
}
//calculate 3 displacement vectors binned by z distance
int findStraindist(struct entry* rd1,struct entry* d1,int size, double * dmir, double deltaz, char * app,int dir)
{
	int count1=0;  
	double * accums; 
	int * counts; 
	int len; 
	double vecs[3]; 
	FILE * fp; 
	char mat[50]; 
	join(mat,"straindist",app,".stat");
	fp=fopen(mat,"w");
	len=(int)(dmir[dir]/deltaz+1); 
	printf("number of bins is %d \n", len); 
	printf("size of bins is %f \n",(float)deltaz); 
	fprintf(fp,"%d \n",len); 
	double zmin; 
	minimum(d1,size,dir,&zmin); //need to determine minimum value of z component
	printf("zmin is %f \n",(float)zmin); 
	counts=(int *)calloc(len*3,sizeof(int));
	accums=(double *)calloc(len*3,sizeof(double));
	if(counts==NULL || accums==NULL )
	{
		printf("no memory \n"); 
		return 0; 
	}
	
	double min=100000;  //determine minimum and maximum displacements to make sure method is working properly 
	double max=-100000;
	printf("binning information \n"); 
	int n, n1, n2, bin; 
	for(n=0; n<size; n++)
	{
		vec(&rd1[n],&d1[n],vecs,dmir); 
		for(n1=0; n1<3; n1++)
		{
			bin=(int)(rd1[n].ent[dir]-zmin)/deltaz; //need to make all z component values positive
			if(bin>=len)
			{
				printf("exceeding size \n");
				continue; //continue going through loop
			}
			counts[bin+n1*len]++;               //increment counter
			accums[bin+n1*len]=accums[bin+len*n1]+vecs[n1]; //accumulate displacements
			if(vecs[n1]>max)
				max=vecs[n1]; 
			if(vecs[n1]<min)
				min=vecs[n1]; 
		}	
	}
	
	printf("outputind binned information \n"); 
	//output binned information
	fprintf(fp,"\n"); 
	for(n=0; n<len*3; n++) //print average of stress tensor at given position
	{
		if(counts[n]!=0)
			accums[n]=accums[n]/counts[n]; 
		count1=counts[n]+count1;
		fprintf(fp,"%f \n",accums[n]);  
	}

	if(count1!=size*3)
	{
		printf("count off, we have %d instead of %d \n",count1,size);
		return 0;  
	}
	fclose(fp);
	printf("minimum displacement is %f and maximum displacement is %f \n",(float)min,(float)max);  
	return 0; 	
}
//we do binning in two dimensions instead of one, so we can produce a 2D color plot/vector plot of stress distribution around steps
//dir1, dir2 tell what directions we are binning
int findStressdistxz(struct entry* d1,int size, double * dmir, double * shear, double deltaz, char * app, int dir1, int dir2)
{

	int count=0;  
	double * accums;
	double deltaz1,deltaz2; 
	deltaz1=deltaz;
	deltaz2=deltaz;  
	int * counts;  
	int len1, len2; 
	FILE * fp; 
	char mat[50]; 
	join(mat,"stressdistxz",app,".stat");
	fp=fopen(mat,"w");
	len1=(int)(dmir[dir1]/deltaz1+1);
	len2=(int)(dmir[dir2]/deltaz2+1);
	printf("number of bins in dimensions 1 and 2 is %d %d \n", len1,len2); 
	printf("size of bins xz is %f \n",(float)deltaz); 
	fprintf(fp,"%d \n",len1);
	fprintf(fp,"%d \n",len2);
	double x1min, x2min; 
	minimum(d1,size,dir1,&x1min); //need to determine minimum value of x1 component
	minimum(d1,size,dir2,&x2min); //need to determine minimum value of x2 component
	printf("x1min is %f, x2min is %f \n",(float)x1min,(float)x2min); 
	counts=(int *)calloc(len1*len2*7,sizeof(int));
	accums=(double *)calloc(len1*len2*7,sizeof(double));
	if(counts==NULL || accums==NULL )
	{
		printf("no memory \n"); 
		return 0; 
	}
	printf("binning information \n"); 
	int n, n1, n2, bin1, bin2; 
	for(n=0; n<size; n++)
	{
	    //remove shear from dimensions so we don't have any problem with non-orthogonal cell
	    	d1[n].ent[2]=d1[n].ent[2];  
		d1[n].ent[1]=d1[n].ent[1]-d1[n].ent[2]*shear[1]/dmir[2]; 
		d1[n].ent[0]=d1[n].ent[0]-d1[n].ent[2]*shear[2]/dmir[2]-d1[n].ent[1]*shear[0]/dmir[1];
	
	    //we bin values for 6 components of stress tensor and potential energy
		for(n1=0; n1<7; n1++)
		{
			if(n1==6)
				n2=PE; 
			else
				n2=SX+n1; 	
			
			
			bin1=(int)(d1[n].ent[dir1]-x1min)/deltaz1; //calculate bin dimensions of given atom  
			bin2=(int)(d1[n].ent[dir2]-x2min)/deltaz2; //calculate bin dimensions of given atom 
			if(bin1>=len1 || bin2>=len2)
			{
				printf("exceeding size \n");
				continue; //continue going through loop, we just don't read in this atom
			} 
			
			
			counts[bin2+bin1*len2+n1*len1*len2]++;               //address scheme: bin1, bin2, data number
			accums[bin2+bin1*len2+n1*len1*len2]=accums[bin2+bin1*len2+n1*len1*len2]+d1[n].ent[n2]; //add stress
		}	
	}
	
	printf("outputind binned information \n"); 
	//output binned information
	fprintf(fp,"\n"); 
	for(n=0; n<len1*len2*7; n++) //print average of stress tensor at given position
	{
		//probably need to update this portion of the code later
		if(counts[n]!=0)
			accums[n]=accums[n]/counts[n]; 
		count=counts[n]+count;
		fprintf(fp,"%f \n",accums[n]);  
	}

	for(n=0; n<len1*len2; n++) //print atom count in each slice
	{
		fprintf(fp,"%d \n",counts[n]);  
	}
	
	fprintf(fp,"%f %f %f \n",dmir[0],dmir[1],dmir[2]); //print xyz dimensions to be used in conversion

	if(count!=size*7)
	{
		printf("count off, we have %d instead of %d \n",count,size);
		return 0;  
	}
	fclose(fp);

        /*
	double shearnew[3];
	shearnew[0]=0; 
	shearnew[1]=0; 
	shearnew[2]=0; 
	 
        writeoutall("test.dump",d1,size,dmir,shearnew,"1"); */

	return 0; 	
}
int findStressdist(struct entry* d1,int size, double * dmir, double deltaz, char * app, int dir)
{
	int count1=0;  
	double * accums; 
	int * counts1; 
	int * counts2; 
	int len; 
	FILE * fp; 
	char mat[50]; 
	join(mat,"stressdist",app,".stat");
	fp=fopen(mat,"w");
	len=(int)(dmir[dir]/deltaz+1); 
	printf("number of bins is %d \n", len); 
	printf("size of bins is %f \n",(float)deltaz); 
	fprintf(fp,"%d \n",len); 
	double zmin; 
	minimum(d1,size,dir,&zmin); //need to determine minimum value of z component
	printf("zmin is %f \n",(float)zmin); 
	counts1=(int *)calloc(len*7,sizeof(int));
	counts2=(int *)calloc(len*7,sizeof(int));
	accums=(double *)calloc(len*7,sizeof(double));
	if(counts1==NULL || counts2==NULL || accums==NULL )
	{
		printf("no memory \n"); 
		return 0; 
	}
	printf("binning information \n"); 
	int n, n1, n2, bin; 
	for(n=0; n<size; n++)
	{
		for(n1=0; n1<7; n1++)
		{
			if(n1==6)
				n2=PE; 
			else
				n2=SX+n1; 	
			
			bin=(int)(d1[n].ent[dir]-zmin)/deltaz; //need to make all z component values positive
			if(bin>=len)
			{
				printf("exceeding size \n");
				continue; //continue going through loop
			} 
			
			if(d1[n].ent[TYPE]==1)
			{
				counts1[bin+n1*len]++;               //increment counter
			}
			else
			{
				counts2[bin+n1*len]++;               //increment counter
			}
			accums[bin+n1*len]=accums[bin+len*n1]+d1[n].ent[n2]; //add stress
		}	
	}
	
	printf("outputind binned information \n"); 
	//output binned information
	fprintf(fp,"\n"); 
	for(n=0; n<len*6; n++) //print average of stress tensor at given position
	{
		//probably need to update this portion of the code later
		if((counts1[n]+counts2[n])!=0)
			accums[n]=accums[n]/(counts1[n]+counts2[n]); 
		count1=counts1[n]+counts2[n]+count1;
		fprintf(fp,"%f \n",accums[n]);  
	}
	//fprintf(fp,"energy \n"); 
	for(n=len*6; n<len*7; n++) //print total energy at given slice
	{
		count1=counts1[n]+counts2[n]+count1;
		fprintf(fp,"%f \n",accums[n]);  
	}
	//fprintf(fp,"bin sizes \n"); 
	//print out Cu atom counts in each slice
	for(n=0; n<len; n++) //print atom count at given slice
	{
		fprintf(fp,"%d \n",counts1[n]);  
	}
	//print out Nb atom counts in each slice
	for(n=0; n<len; n++) //print atom count at given slice
	{
		fprintf(fp,"%d \n",counts2[n]);  
	}
	
	fprintf(fp,"%f %f %f \n",dmir[0],dmir[1],dmir[2]); //print xyz dimensions to be used in conversion

	if(count1!=size*7)
	{
		printf("count off, we have %d instead of %d \n",count1,size);
		return 0;  
	}
	fclose(fp); 
	return 0; 	
}

//need reference lattice and deformed lattice as input arguments, and reference number, output of code is stress and strain fields as function of z
//each of output files can be zipped together or whatever for post processing by matlab
int main(int argc, char* argv[])
{
	char * processfile; 
	char * referencefile;
	char * command;   
	int refnumber, dir, type,isxztest,stepnumber,n,includestressfields;  
	double deltaz,z1,z2,stepheight; 
	double steplocations[20]; 
	
	referencefile=argv[1];  //reference file
	processfile=argv[2]; //file to be processes 
	refnumber=(int)atof(argv[3]);  //number to be appended to output files
	
	command=argv[4]; 

	//get additional input arguments depending on which command was specified
	if(strcmp(command,"stressdist")==0)
	{
		deltaz=(double)atof(argv[5]); //get size of slices
		printf("stress dist command with %f for slices \n",(float)deltaz); 
	}	
	if(strcmp(command,"stressdistxz")==0)
	{
		deltaz=(double)atof(argv[5]); //get size of slices
		printf("stressdistxz command with %f for slices \n",(float)deltaz); 
	}	
	if(strcmp(command,"straindist")==0)
	{
		deltaz=(double)atof(argv[5]); //get size of slices
		printf("strain dist command with %f for slices \n",(float)deltaz); 
	}
	if(strcmp(command,"stressdistdir")==0)
	{
		deltaz=(double)atof(argv[5]); //get size of slices
		printf("stress dist command with %f for slices \n",(float)deltaz);
		dir=(int)atof(argv[6]); 
	}
	if(strcmp(command,"straindistdir")==0)
	{
		deltaz=(double)atof(argv[5]); //get size of slices
		printf("strain dist command with %f for slices \n",(float)deltaz);
		dir=(int)atof(argv[6]);  
	}
	if(strcmp(command,"minmax")==0)
	{
		printf("calculating minimum maximum components \n"); 
	}
	if(strcmp(command,"writeout")==0)
	{
		printf("reading in and then printing out dumpfile \n"); 
	}
	if(strcmp(command,"reffile")==0)
	{
		z1=(double)atof(argv[5]); 
		z2=(double)atof(argv[6]); 
		type=(int)atof(argv[7]);
		if(type==4)  //we need to get several additional pieces of information about step location, etc. 
		{
			stepnumber=(int)atof(argv[8]);
			printf("%d stepnumber \n",stepnumber); 
			for(n=0; n<(stepnumber*2); n++)
			{
				steplocations[n]=(double)atof(argv[9+n]); 
				printf("%f \n", steplocations[n]); 
			}
			stepheight=(double)atof(argv[9+stepnumber*2]);
			printf("%f \n",stepheight);  
		} 
		printf("calculating reference file for disregistry calculation z range %lf %lf type %d \n",z1,z2,type); 
	}
	if(strcmp(command,"avgz")==0)
	{
		printf("calculating average z direction lattice spacing \n"); 
	}
	if(strcmp(command,"disregistry")==0)
	{
		type=(int)atof(argv[5]); 
		includestressfields=(int)atof(argv[6]); //determine whether including stress fields of Nb layer at interface
		printf("calculating disregistry type %d \n",type); 
	}
	
	struct entry * ref; //reference atoms
	struct entry * data; //data atoms	
		
	int atoms; 
	double dims[3]; 
	double shear[3]; 

	//set append word for files
	char append [10]; 
	sprintf(append,"%d",refnumber); //for now just have append reference number
	printf("append is %s \n",append); 	

	printf("Reading in reference atoms \n"); 
	n=readinall(referencefile,&ref,&atoms,dims,shear);
	if(n==-1)
	{
		printf("exiting program, bad readin all \n"); 
		return 0; 
	} 
	printf("Reading in data atoms \n"); 
	n=readinall(processfile,&data,&atoms,dims,shear);
	if(n==-1)
	{
		printf("exiting program, bad readin all \n"); 
		return 0; 
	} 

	
	//now do stuff with the data based on the input command
 	if(strcmp(command,"stressdist")==0)
	{
		findStressdist(data,atoms,dims,deltaz,append,2); //as default uses z component
	}
 	if(strcmp(command,"straindist")==0)
	{
		findStraindist(ref,data,atoms,dims,deltaz,append,2); //as default uses z component
	}
 	if(strcmp(command,"stressdistdir")==0)
	{
		findStressdist(data,atoms,dims,deltaz,append,dir); //as default uses z component
	}
	if(strcmp(command,"stressdistxz")==0)
	{
		findStressdistxz(data,atoms,dims,shear,deltaz,append,1,2); //we want to bin-data by x z components so we can visualize 2D field 
	}
 	if(strcmp(command,"straindistdir")==0)
	{
		findStraindist(ref,data,atoms,dims,deltaz,append,dir); //as default uses z component
	}
	if(strcmp(command,"avgz")==0)
	{
		calculateAvgZ(data,atoms,dims); 
	}
	if(strcmp(command,"minmax")==0)
	{
		printf("calculating minimum maximum components \n");
		double min; 
		double max; 
		int n; 
		for(n=0; n<3; n++)
		{
			minimum(ref,atoms,n,&min); 
			maximum(ref,atoms,n,&max);
			printf("min and max for %d component is %f %f \n",n,(float)min,(float)max);  
		}
 
	}
	if(strcmp(command,"writeout")==0)
	{
		writeoutall(referencefile,ref,atoms,dims,shear,append); 
	}
	if(strcmp(command,"reffile")==0)
	{
		createRefFile(ref,atoms,dims,shear,z1,z2,type,stepnumber,steplocations,stepheight,referencefile,append); 
	}
	if(strcmp(command,"disregistry")==0)
	{
		calculateDisregistry(ref,atoms,dims,shear,type,includestressfields,append); 
	}


	free(ref); //free allocations for reference and data atoms
	free(data); 
}
