
#include <stdio.h>
#include <math.h>

int main()
{
	//FILE * fp=fopen("CuNb.KS1min.data","r"); 
	FILE * fp=fopen("CuNb_KS1.data","r"); 
	
	if(fp==NULL)
	{
		printf("invalid file \n");
		return;
	}

	int n; 	

	char bob [100]; 

	for(n=0; n<19; n++)
		fgets(bob,100,fp); //ignore input lines

	int count1,count2,check,num1,num2; 
	float pos; 

	count1=0; 
	count2=0;
 
	while(1)
	{
		check=fscanf(fp,"%d %d %f %f %f \n",&num1,&num2,&pos,&pos,&pos); 
		if(check==EOF)
		{
			break;
		}
		else
		{
			if(num2==1)
				count1++; 
			else
				count2++; 	
		} 
	}


	printf("have %d cu atoms and %d nb atoms for a total of %d \n",count1,count2,count1+count2); 

	fclose(fp); 
	return 0; 
}
