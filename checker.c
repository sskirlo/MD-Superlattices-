
//program determines which ids are missing from data file

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
	FILE * fp=fopen("CuNb_KS1min.data","r"); 
	FILE * missing=fopen("missing.txt","w"); 

	if(fp==NULL)
	{
		printf("invalid file \n");
		return;
	}

	int len=57144; 

	int * bob; 
	bob=(int *)calloc(len,sizeof(int));
	if(bob==NULL)
	{
		printf("bad memory call \n"); 
	}
	int n;  
	for(n=0; n<len; n++)
	{
		bob[n]=-1; 
	}

	int type,num;
	float pos;

	char text [100];
	for(n=0; n<19; n++)
	{
		fgets(text,100,fp);
		printf("%s ",text);  
	}
	for(n=0; n<len; n++)
	{
		fscanf(fp,"%d %d %f %f %f ",&num,&type,&pos,&pos,&pos);
		//printf("num %d \n",num); 
		bob[num-1]=1;  	
	}

	for(n=0; n<len; n++)
	{
		if(bob[n]==-1)
		{
			fprintf(missing,"%d \n",n+1); 
		}
	}

	fclose(missing); 
	fclose(fp); 
	return 0; 
}
