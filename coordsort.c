
#include <stdio.h>
#include <math.h>

struct bob
{
	int id; 
	int coordination; 
};
//function compares based on coordination number, sorted this way after entering
int comparecoords(const void *ptr1, const void *ptr2) //compare function to be passed to qsort
{
	int c1, c2; 
	struct bob* p1;
	struct bob* p2; 
	p1=(struct bob*)ptr1; 
	p2=(struct bob*)ptr2;  
	c1=(*p1).coordination; 
	c2=(*p2).coordination;
	if(c1<c2)
		return 1;
	else if(c1>c2)
		return -1;
	else
		return 0;
}
int main()
{
	FILE * fp=fopen("cunb_view0.dump","r"); 
	FILE * fout=fopen("vacids.txt","w"); 	

	if(fp==NULL)
	{
		printf("invalid file \n");
		return;
	}

	int n; 	

	char bob [100];

	//scan in number of atoms
	int num; 
	fscanf(fp,"%s %s %s %s %d ",&bob,&bob,&bob,&bob,&num);  

	for(n=0; n<17; n++)
		fgets(bob,100,fp); //ignore input lines

	
	float x; 

	struct bob * joe; 
	joe=(struct bob *)calloc(num,sizeof(struct bob)); 
		

	for(n=0; n<num; n++)
	{
		fscanf(fp,"%f %f %f %d %d \n",&x,&x,&x,&(joe[n].id),&(joe[n].coordination)); 
	}
	qsort((void*)joe,num,sizeof(struct bob),comparecoords); 
	
	for(n=0; n<num; n++)
	{
		//fprintf(fout,"%d %d \n",joe[n].id,joe[n].coordination); 	
		fprintf(fout,"%d \n",joe[n].id); 	
	}
		
	fclose(fp); 
	fclose(fout); 
	return 0; 
}
