#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<omp.h>

#define FLN 1000

void main(int argv, char *argc[])
{
	int first,last;

	first = atoi(argc[1]);
	last = atoi(argc[2]);

	int i,id,sym;
	char buff[1000];

#pragma omp parallel
{
#pragma omp for private(buff)     // we need to use private, otherwise the buff is shared
	for(i=first;i<=last;i++)
	{
		sprintf(buff,"./ChemNetworks.exe Input water%d.xyz \n",i);  // this is the command of running ChemNetworks
		sym=system(buff);   // this is how C code can call a command
		id=omp_get_thread_num();
		printf("  %d from proc %d return %d\n",i,id,sym);
	}
}

}


