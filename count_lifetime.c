#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define FLN 100

int main(int argc , char *argv[])
{
	FILE *fp,*fpnew;
	long i=0;
	long j=0;
	long max=1;
	long totnum,totlftm,totsquare;
	double avglftm,stddev;
	char buffer[10];

	char finput[FLN];      /* intput file name */
        char foutput[FLN];

// Start reading Input Files
	if(argc != 2){
	printf("Usuage: %s Lifetimes.DIMER\n", argv[0]);
	exit(-1);
	}

	sprintf(finput, "%s", argv[1]);

	if((fp = fopen(finput,"r")) == NULL){
		printf("Cannot open file %s\n", finput);
		exit(-1);
	}

	sprintf(foutput, "lifetime-count");
	fpnew=fopen(foutput, "w");	

	while(fscanf(fp,"%ld",&i)!=EOF)
	{
		fgets(buffer,sizeof(buffer),fp);
		if(i>max)
			max=i;
	}
	printf("-----max=:%ld\n", max);

	long NUM[max+1];
	long TNUM;

	for(i=0;i<max+1;i++)
		NUM[i]=0;
	TNUM=0;

	printf("counting durations\n");
	rewind(fp);
	while(fscanf(fp, "%ld", &i)!=EOF)
	{
		if(fgets(buffer,sizeof(buffer),fp)==NULL)
			break;
		if(i>0 && i <=max)
		{
			NUM[i]++;
			TNUM++;
		}
	}

	totnum=totlftm=totsquare=0;
	avglftm=0.0;

	printf("exporting\n");
	for(i=1;i<max+1;i++)
	{
		if(NUM[i]!=0)
		{
			totnum += NUM[i];
			totlftm += i*NUM[i];
			totsquare += i*i*NUM[i];
			fprintf(fpnew, "duration %ld number %ld percent %f\n", i,NUM[i],(NUM[i]+0.0)/TNUM);
		}
	}
	avglftm = (totlftm+0.0)/totnum;
	printf("%d %d %d\n",totnum,totlftm,totsquare);
	stddev = sqrt( (totsquare+0.0)/(totnum+0.0) - (avglftm+0.0)*(avglftm+0.0) );
	printf("\naverage lifetime of %s is: %.5f , standard-deviation %.5f (snaps)\n\n", argv[1], avglftm, stddev);
	fprintf(fpnew,"\n\n\naverage lifetime of %s is: %.5f snaps, standard deviation: %.5f snaps\n\n", argv[1], avglftm, stddev);

	fclose(fp);
	fclose(fpnew);
}


