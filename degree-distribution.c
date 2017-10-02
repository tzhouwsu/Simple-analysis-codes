/* 
by Tiecheng Zhou: tiecheng.zhou@email.wsu.edu -- 2014.0127

this code can:
1)calculate degree distribution averaged over snapshots and water-ids

need all graph files: GraphGeod files
*/


#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define Nnodes 216
#define NSNAP 40000	//total number of snapshot

int degreecount(long a[Nnodes], int snap)
{
	long deg;
	char filename[100];
	FILE *fpp;
	char temp[100];
	int i,j;

	sprintf(filename, "Water.input.water%d.xyz.water%d.xyz.GraphGeod",snap,snap);
	if((fpp=fopen(filename, "r"))==NULL)
	{
		printf("ERROR: can't find the graph files for degree analysis\n");
		return(-1);
	}
	
	int counted[Nnodes][Nnodes]; // matrix indicating a particular H-bond i-j is counted within the same graph
	for(i=0;i<Nnodes;i++)
		for(j=0;j<Nnodes;j++)
			counted[i][j]=0;

	deg=0;
	while(fscanf(fpp, "%d %d ", &i, &j)!=EOF)
	{
		if(fgets(temp, sizeof(temp), fpp)==NULL)
			break;
		if(counted[i-1][j-1]==0) // H-bond i-j has not been counted
		{
			a[i]++;
			a[j]++;
			counted[i-1][j-1]=counted[j-1][i-1]=1; // set to 1 to avoid double counting
		}
	}

	fclose(fpp);
	return(0);
}

int main()
{
	long degree[10];	//from degree 0 to 9
	long snap[Nnodes+1];	//degree of each node at a given snapshot
	int m,n,t;
	int sym=0;
	float percent[10];

//initialize total degree numbers
	for(m=0;m<10;m++)
		degree[m]=0;

	for(t=1;t<=NSNAP;t++)
	{
		//initilize degree of nodes
		for(n=1;n<=Nnodes;n++)
			snap[n]=0;
		
		printf("counting snap %d\n", t);
	
		sym=degreecount(snap,t);
		if(sym==-1)
			break;

		for(n=1;n<=Nnodes;n++)
		{
			degree[snap[n]]+=1;	//node n has degree of snap[n], so degree[snap[n]]+1
		}
	}
	if(sym==-1)
	{
		printf("\n----ERROR-----\n");
		return(-1);
	}

	for(m=0;m<10;m++)
		percent[m]=(degree[m]+0.0)/(Nnodes*NSNAP);

	printf("\ntotal number of degree distribution:\n");
	printf("-------------------------------------------\n");
	for(m=0;m<10;m++)
	{
		printf("degree %d number %lu Percent = %f\n", m, degree[m],100*percent[m]);
	}
	
	return(0);

}

