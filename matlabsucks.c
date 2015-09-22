#include <stdio.h>
#include <stdlib.h>

//order of arguments is input filename, output filename, and readin count
int main(int argc,char * argv[])
{
	FILE *fpin; 
	FILE *fpout; 
	fpin=fopen(argv[1],"r"); 
	fpout=fopen(argv[2],"w"); 

	if(fpin==NULL)
	{
		printf("bad input file name \n"); 
		return 0;
	} 
	if(fpout==NULL)
	{
		printf("bad output file name \n"); 
		return 0;
	} 

	int counts; 
	counts=(int)atof(argv[3]); 	

	int buff=300; 
	char text [300]; 

	int n,n2; 
	float d1,d2; 
	for(n=0; n<counts; n++)
	{
		for(n2=0; n2<5; n2++)
		{			
			fgets(text,buff,fpin); 
		}	
		
		//fgets(text,buff,fpin); 
		//printf("%s \n",text); 
		fscanf(fpin,"%s %f %f ",text,&d1,&d2); 
		fprintf(fpout,"%f \n",d2); 
		fgets(text,buff,fpin); //need to get end line
		fgets(text,buff,fpin); 
	}

	fclose(fpin); 
	fclose(fpout); 
	return 0; 
}
