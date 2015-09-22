
//determines minimum and maximum coordinates in data file

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
	FILE * fp=fopen("CuNb_KS1.data","r"); 

	if(fp==NULL)
	{
		printf("invalid file \n");
		return;
	}

	int len=57144; 

	int type,num;
	float pos[3];

	char text [100];
	for(n=0; n<19; n++)
	{
		fgets(text,100,fp);
		printf("%s ",text);  
	}
	float min=100000; 
	float max=-100000; 

	for(n=0; n<len; n++)
	{
		fscanf(fp,"%d %d %f %f %f ",&num,&type,&(pos[0]),&(pos[1]),&(pos[2]));
		if(pos[2]>max)
			max=pos[2]; 
		if(pos[2]<min)
			min=pos[2]; 
	}

	printf("max is %f, min is %f \n",max,min); 

	fclose(fp); 
	return 0; 
}
