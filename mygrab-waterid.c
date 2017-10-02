
// this code is to grab the water coordinates based on several molecule id (user entered)
//    need the xyz coordinates of all the waters, and I should know a list of water_ids that need to be grabed
//    also consider PBC (need PBC box size), shift all the coordinate so that the grabed coordinates are in the same PBC
//  by Tiecheng Zhou

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define FLN 1000
#define Buffsize 10000

int pbcshift(double result[3],double xi,double yi,double zi,double refx,double refy,double refz, double Xsize,double Ysize,double Zsize)
{
	result[0]=result[1]=result[2]=-1.0;
	if(xi-refx > Xsize/2)
		result[0]=xi-Xsize;
	else if(xi-refx < -1.0*Xsize/2)
		result[0]=xi+Xsize;
	else
		result[0]=xi;

	if(yi-refy > Ysize/2)
		result[1]=yi-Ysize;
	else if(yi-refy < -1.0*Ysize/2)
		result[1]=yi+Ysize;
	else
		result[1]=yi;

	if(zi-refz > Zsize/2)
		result[2]=zi-Zsize;
	else if(zi-refz < -1.0*Zsize/2)
		result[2]=zi+Zsize;
	else
		result[2]=zi;

//printf(" %f %f %f %f %f %f %f %f %f\n",result[0],result[1],result[2],xi,yi,zi,refx,refy,refz);

	return(0);
}

double getdistance(double xi,double yi,double zi,double refx,double refy,double refz, double Xsize,double Ysize,double Zsize)  // calculate the distance of two atoms
{
	double distx,disty,distz,result;
	if(fabs(xi-refx) > Xsize/2)
		distx = Xsize - fabs(xi-refx);
	else
		distx = fabs(xi-refx);
	
	if(fabs(yi-refy) > Ysize/2)
		disty = Ysize - fabs(yi-refy);
	else
		disty = fabs(yi-refy);
	
	if(fabs(zi-refz) > Zsize/2)
		distz = Zsize - fabs(zi-refz);
	else
		distz = fabs(zi-refz);
	
	result = sqrt(distx*distx + disty*disty + distz*distz);

//	printf(" %f %f %f %f %f %f : %f\n",xi,yi,zi,refx,refy,refz,result);

	return(result);
}


int main(int argc, char *argv[])
{
	// here is pre-run input
	char filename1[FLN],filename2[FLN];	
	FILE *fip1,*fip2;
	double sizeX,sizeY,sizeZ;
	int key,dum;
	int Nwaters,nindex;
	int *list_index;

	printf("Enter the input coordinates that has all water molecules:\n");
	scanf("%s",filename1);
	if( (fip1=fopen(filename1,"r"))==NULL )
	{
		printf("  Error: cannot find the input file %s\n",filename1);
		exit(-1);
	}
	if( fscanf(fip1,"%d",&Nwaters) != 1)
	{
		printf(" Error: did not read 'Nwaters' from file %s\n",filename1);
		fclose(fip1);
		exit(-1);
	}
	if(Nwaters <= 0)
	{
		printf(" wrong total number of waters %d in file %s\n",Nwaters,filename1);
		fclose(fip1);
		exit(-1);
	}
	else
	{
		Nwaters = Nwaters/3;
		list_index = (int *) calloc(Nwaters+1,sizeof(int));
	}

	printf("Enter the box size:\n  Xsize = (use floating number)\n");
	scanf("%lf",&sizeX);
	printf("  Xsize = %lf \n",sizeX);
	printf("Enter the box size:\n  Ysize = (use floating number)\n");
	scanf("%lf",&sizeY);
	printf("  Ysize = %lf \n",sizeY);
	printf("Enter the box size:\n  Zsize = (use floating number)\n");
	scanf("%lf",&sizeZ);
	printf("  Zsize = %lf \n",sizeZ);

	printf("Enter integer 1 if you have a file for the water index to read\n");
	printf("         or   0 if you want to enter by input (end by 0)\n");
	scanf("%d",&key);

	if(key == 1)
	{
		printf("  now enter the file-name (contains integer numbers)\n");	
		scanf("%s",filename2);
		if( (fip2=fopen(filename2,"r"))==NULL )
		{
			printf("  error: cannot find the file %s\n",filename2);
			fclose(fip1);
			exit(-2);
		}

	}
	else if(key == 0)
	{
		printf("  now enter the water_index(s) you want to grab: (end by entering 0)\n");
		dum=1; nindex=0;
		while(dum != 0)
		{
			scanf("%d",&dum);
			if(dum >= Nwaters)
			{
				printf("  wrong number, %d larger than total number %d\n",dum,Nwaters);
				break;
			}
			else if(dum > 0)
			{
				nindex++;
				list_index[nindex]=dum;
			}
		}
		if(nindex == 0)
		{
			fclose(fip1);
			exit(-3);
		}
	}
	else
	{
		printf(" error: did not get integer 1/0, skip the code now\n");
		fclose(fip1);
		exit(-3);
	}

	// start to read the all-water coordinate file
	double Crdw[Nwaters+1][9];  // the coordinates of waters
	char buffer[Buffsize];
	int num,i,j;
	char lab;
	double x,y,z;
	for(i=0;i<=Nwaters;i++)
		for(j=0;j<9;j++)
			Crdw[i][j]=0.0;

	fgets(buffer,sizeof(buffer),fip1); 
	fgets(buffer,sizeof(buffer),fip1);
	num=0;
	while(fscanf(fip1,"%c %le %le %le",&lab,&x,&y,&z)==4 && num<Nwaters)
	{
		fgets(buffer,sizeof(buffer),fip1);
		if(lab=='O')   // the water.xyz should be in a specific order that Oxygen coordinates is listed prior to Hydrogens
		{
			num++;
			Crdw[num][0] = x;
			Crdw[num][1] = y;
			Crdw[num][2] = z;

			fscanf(fip1,"%c %le %le %le",&lab,&x,&y,&z);
			fgets(buffer,sizeof(buffer),fip1);
			if(lab=='H')
			{
				Crdw[num][3] = x;
				Crdw[num][4] = y;
				Crdw[num][5] = z;
			}
			else
				break;

			fscanf(fip1,"%c %le %le %le",&lab,&x,&y,&z);
			fgets(buffer,sizeof(buffer),fip1);
			if(lab=='H')
			{
				Crdw[num][6] = x;
				Crdw[num][7] = y;
				Crdw[num][8] = z;
			}
			else
				break;

		}
		else if(lab=='H')  
			break;

//		printf(" %d %f %f %f %f %f %f %f %f %f\n",num,Crdw[num][0],Crdw[num][1],Crdw[num][2],Crdw[num][3],Crdw[num][4],Crdw[num][5],Crdw[num][6],Crdw[num][7],Crdw[num][8]);
	}

	if(num != Nwaters)
	{
		printf("Error: in reading the %s : %d %d\n",filename1,num,Nwaters);
		fclose(fip1);
		exit(-3);
	}
	else
		fclose(fip1);


	// start to read the water_index from a file
	if(key == 1)
	{
		num=0;
		while(fscanf(fip2,"%d",&dum)==1 && num<Nwaters)
		{
//	printf("%d %d\n",num,dum);
			if(dum < 0 || dum > Nwaters)
			{
				fclose(fip2);
				exit(-4);
			}
			else
			{
				num++;
				list_index[num]=dum;
			}
		}

		nindex = num;
		fclose(fip2);
	}

	// start to grab the coordinate
	FILE *fop;
	double refx,refy,refz;
	double Ox,Oy,Oz,H1x,H1y,H1z,H2x,H2y,H2z;
	double shiftresult[3];

	fop = fopen("output.xyz","w");

	fprintf(fop,"   %d  \n\n",nindex*3);
	refx = Crdw[1][0]; refy = Crdw[1][1]; refz = Crdw[1][2];  // take water_1 as the initial reference point (PBC) 
	for(num=1;num<=nindex;num++)
	{
		dum = list_index[num]; // this is the water_index to be grabbed
		Ox = Crdw[dum][0]; Oy = Crdw[dum][1]; Oz = Crdw[dum][2];
		H1x = Crdw[dum][3]; H1y = Crdw[dum][4]; H1z = Crdw[dum][5];
		H2x = Crdw[dum][6]; H2y = Crdw[dum][7]; H2z = Crdw[dum][8];

		pbcshift(shiftresult,Ox,Oy,Oz,refx,refy,refz,sizeX,sizeY,sizeZ);  // shift the oxygen coordinates
		Ox = shiftresult[0]; Oy = shiftresult[1]; Oz = shiftresult[2];
		pbcshift(shiftresult,H1x,H1y,H1z,refx,refy,refz,sizeX,sizeY,sizeZ);  
		H1x = shiftresult[0]; H1y = shiftresult[1]; H1z = shiftresult[2];
		pbcshift(shiftresult,H2x,H2y,H2z,refx,refy,refz,sizeX,sizeY,sizeZ);  
		H2x = shiftresult[0]; H2y = shiftresult[1]; H2z = shiftresult[2];

		fprintf(fop,"O  %lf  %lf  %lf\n",Ox,Oy,Oz);
		fprintf(fop,"H  %lf  %lf  %lf\n",H1x,H1y,H1z);
		fprintf(fop,"H  %lf  %lf  %lf\n",H2x,H2y,H2z);

		refx = Ox; refy = Oy; refz = Oz;  // reset the reference to the current water, (if you entered a water HB chain, then every connected water should be in the same PBC
	}

	fclose(fop);

	return(0);
	
}





