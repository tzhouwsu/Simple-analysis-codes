
/* this code is to get the statistic histogram for data in the inputfile
 * the analyzed data is in column x in the inputfile (x need to be specified)
 * also the lower-boundary, the upper-boundary, and the binsize need to be specified
 */


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define FLN 1000000

int main(int argc, char *argv[])
{
	char filename[FLN];
	float istart,iend,ibin;
	float min,max;
	int icolumn,nbins,nskip;
	FILE *fip;

	if(argc != 7 )
	{
		printf("How to use it:\n  %s input-file n_skip column start end bin-size\n\n",argv[0]);
		exit(-1);
	}
	sprintf(filename,"%s",argv[1]);
	if((fip=fopen(filename,"r"))==NULL)
	{
		printf("Error: Can not find the input-file\n");
		exit(-1);
	}
	nskip = atoi(argv[2]);
	icolumn = atoi(argv[3]);
	istart = atof(argv[4]);
	iend = atof(argv[5]);
	ibin = atof(argv[6]);
	nbins = (int) ((iend-istart)/ibin);

	printf("Get statistics for file %s\n with column %d and range %f %f %f (skip first %d lines)\n\n",argv[1],icolumn,istart,iend,ibin,nskip);

	int *Obsnum;
	int i,j,line,id;
	double value,total;
	char buff[FLN],*token;

	Obsnum = (int *)malloc((nbins+1)*sizeof(int));

	line=0;total=0.0;
	for(i=1;i<=nskip;i++)
		fgets(buff,sizeof(buff),fip);
	while(fgets(buff,sizeof(buff),fip) != NULL)   // read line by line
	{
		line++;
//printf("%d %s\n",line,buff);
		token = strtok(buff," \n");  // the columns are separated by space ' '
		for(i=2;i <= icolumn;i++)
			token = strtok(NULL," \n");

		value = atof(token);   // here I assume the value is in floating number
		id = (int) ((value-istart)/ibin);
		if(id<0)
			id=0;
		if(id>nbins)
			id=nbins;

		Obsnum[id] += 1;
		total += value;

		if(line ==  1)
			min=max=value;
		else
		{
			if(value < min)
				min = value;
			if(value > max)
				max = value;
		}

	}

	printf("Total lines %d - Print the output\n\n",line);
	for(i=0;i<nbins;i++)
	{
		printf("%f %d\n",istart+i*ibin,Obsnum[i]);
	}

	printf("\naverage value for %d lines: %f\n",line,total/(line+0.0));
	printf("  min value %f , max value %f\n\n",min,max);

	fclose(fip);
	return(0);

}


