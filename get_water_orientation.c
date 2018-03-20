/* this is to calculate the water dipole orientation with respect to interfacial normal at water:hexane interface */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

#define Boxx 40.1755221
#define Boxy 40.179491
#define Boxz 141.635303
#define Natoms 21335
#define Nsnaps 20000
#define Nbins 1400     // change it for diff. resolution in z
#define PI 3.14159265 


float mycaltheta(float Ox, float Oy, float Oz, float H1x, float H1y, float H1z, float H2x, float H2y, float H2z)
{
	float result=0.0;
	float dipx,dipy,dipz,length;

	if(H1x-Ox > Boxx/2)
		H1x = H1x - Boxx;
	else if(H1x-Ox < -1*Boxx/2)
		H1x = H1x + Boxx;

	if(H1y-Oy > Boxy/2)
		H1y = H1y - Boxy;
	else if(H1y-Oy < -1*Boxy/2)
		H1y = H1y + Boxy;

	if(H1z-Oz > Boxz/2)
		H1z = H1z - Boxz;
	else if(H1z-Oz < -1*Boxz/2)
		H1z = H1z + Boxz; 

	if(H2x-Ox > Boxx/2)
		H2x = H2x - Boxx;
	else if(H2x-Ox < -1*Boxx/2)
		H2x = H2x + Boxx;

	if(H2y-Oy > Boxy/2)
		H2y = H2y - Boxy;
	else if(H2y-Oy < -1*Boxy/2)
		H2y = H2y + Boxy;

	if(H2z-Oz > Boxz/2)
		H2z = H2z - Boxz;
	else if(H2z-Oz < -1*Boxz/2)
		H2z = H2z + Boxz; 

	dipx = H1x+H2x-Ox*2;
	dipy = H1y+H2y-Oy*2;
	dipz = H1z+H2z-Oz*2;
	length = sqrt(dipx*dipx+dipy*dipy+dipz*dipz);
	result = dipz/length;  // the normal of the interface is along z direction
//	printf("%f %f %f %f %f %f %f %f %f %f %f %f %f\n",Ox,Oy,Oz,H1x,H1y,H1z,H2x,H2y,H2z,dipx,dipy,dipz,result);

	return(result);
}

int mycheckinterfacial(int snapid, int moleculeid)  // moleculeid is the molecule id of water from trajectory file
{
	int idconvert; // the moleculeid from the trajectory file is different from that in the GraphGeod file
	int result;
/*
note-2016-07-11
the molecule id from trajectory file is different from that from GraphGeod file
 from trajectory
    Number of waters: 3205
        id: 294 --> 3498
    Number of hexanes: 586
        id: 1 --> 293, 3499 --> 3791
 from GraphGeod
    Number of waters: 3205
        id: 1 --> 3205
    Number of hexanes: 586
        id: 3206 --> 3791
I checked the corrdinates, and they are in increasing order, so just need to shift the molecule id
*/

	idconvert = moleculeid - 293; // idconvert is the molecule id from GraphGeod file

	FILE *fgg;
	char filename[100];
	int id1,id2;

	sprintf(filename,"Input.water%d.xyz.solB%d.xyz.GraphGeod",snapid,snapid);
	if((fgg=fopen(filename,"r"))==NULL)
	{
		printf("Error: in function 'mycheckinterfacial', can not find the GraphGeod file at snap %d\n",snapid);
		result=-1;
	}
	else
	{
		result=0;  // by default, this water molecule (moleculeid,or idconvert) is not interfacial water
		rewind(fgg);
		while(fscanf(fgg,"%d %d",&id1,&id2)==2)
		{
			if(id1==idconvert && (id2>=3206 && id2<=3791))   // there is an interaction between this water and one hexane
				result=1;
			if(id2==idconvert && (id1>=3206 && id1<=3791))
				result=1;
		}

		fclose(fgg);		
	}

	return(result);
}


int main()
{
	int snap;
	FILE *fip,*fp;
	char buff[200];
	int i,j,index,idmol,idatm,ibin,itf;
	float chg,x,y,z,Ox,Oy,Oz,H1x,H1y,H1z,H2x,H2y,H2z;
	float costheta,theta,sz;
	float statang[Nbins+1],countnum[Nbins+1],statcosang[Nbins+1],statsz[Nbins+1];
	float iang[Nbins+1],icount[Nbins+1],icosang[Nbins+1],isz[Nbins+1],bang[Nbins+1],bcount[Nbins+1],bcosang[Nbins+1],bsz[Nbins+1]; // separate to interfacial water and, bulk-like water, 2016.0711
	int nsample,r,k; // sample 20000 waters in bulk-water region, to test whether the dipole fluctuation is caused by small number of observations

	srand(time(NULL));
	k=0;nsample=20000;

	for(i=0;i<=Nbins;i++)
	{
		statang[i]=countnum[i]=statcosang[i]=statsz[i]=0.0;
		iang[i]=icount[i]=icosang[i]=isz[i]=0.0;
		bang[i]=bcount[i]=bcosang[i]=bsz[i]=0.0;
	}

	fp=fopen("result-orient","w");
	if((fip=fopen("/home/alexm/hexane/WATRER_HEXANE/wat-hex-298-nve.lammpstrj","r"))==NULL)
	{
		printf("Can not find the input file\n");
		fclose(fip);
		fclose(fp);
		return(-1);
	}

	rewind(fip);
	for(snap=1;snap<=Nsnaps;snap++)
	{
		printf("for snap %d ...\n",snap);
		for(i=1;i<=9;i++)    // skip the head part
			fgets(buff,sizeof(buff),fip);

		for(i=1;i<=Natoms+1;)
		{
			if(fscanf(fip,"%d %d %d %f %f %f %f",&index,&idmol,&idatm,&chg,&x,&y,&z)!=7)
				break;

			fgets(buff,sizeof(buff),fip); // skip the '\n'
			i++;
			j=0;
			if(idatm==4 && i<=Natoms+1)    // this is an Oxygen atom, then there are two Hydrogens following it
			{
				Ox=x;Oy=y;Oz=z;
				j=1;
				if((fscanf(fip,"%d %d %d %f %f %f %f",&index,&idmol,&idatm,&chg,&x,&y,&z))==7)
				{
					fgets(buff,sizeof(buff),fip);
					i++;
					if(idatm==5 && i<=Natoms+1)   // this is a Hydrogen
					{
						H1x=x;H1y=y;H1z=z;
						j++;
					}
				}
				if((fscanf(fip,"%d %d %d %f %f %f %f",&index,&idmol,&idatm,&chg,&x,&y,&z))==7)
				{
					fgets(buff,sizeof(buff),fip);
					i++;
					if(idatm==5 && i<=Natoms+1)  // this is a Hydrogen
					{
						H2x=x;H2y=y;H2z=z;
						j++;
					}
				}

				if(j==3)  // successfully read the coordinate of a water molecule
				{
					costheta=mycaltheta(Ox,Oy,Oz,H1x,H1y,H1z,H2x,H2y,H2z);
					theta = acos(costheta)*180.0/PI;
					sz = 0.5*(3*costheta*costheta-1);

//	if(Oz >= 57.0 && Oz<=62.0) // this is the layer in the bulk-water limit, 2016.07.18
//	{
//		r=rand()%100; // create a random number between 0 to 99
//		if(r >=50 && k<nsample) // randomly select 1000 waters
//		{
					ibin = (int) (Oz*10);       // change with respect to the z resolution
					statang[ibin] += theta;
					statcosang[ibin] += costheta;
					statsz[ibin] += sz; 
					countnum[ibin] += 1;

					itf = mycheckinterfacial(snap,idmol); // check whether this water molecule is interfacial water or not
					if(itf == 1) // this is an interfacial water
					{
						iang[ibin] += theta;
						icosang[ibin] += costheta;
						isz[ibin] += sz;
						icount[ibin] += 1;
					}
					else if(itf == 0)   // this is a bulk-like water
					{
						bang[ibin] += theta;
						bcosang[ibin] += costheta;
						bsz[ibin] += sz;
						bcount[ibin] += 1;
					}
//			k++;
//		}
//	}

				//	printf("%d %d %f %f %d\n",snap,idmol,Oz,costheta,ibin);
				}
			}
			
		}

	}

	printf("get statistics ..\n");
	fprintf(fp,"\n---for all the waters---\n\n");
	fprintf(fp,"z avg-ang avg-cosang avg-sz\n\n");
	for(i=0;i<=Nbins;i++)
	{
		if(countnum[i] == 0)
			fprintf(fp,"%.3f %f %f %f %.1f\n",i/10.0+0.05,0.0,0.0,0.0,0.0);   // change for diff. z resolution
		else if(countnum[i] > 0)
			fprintf(fp,"%.3f %f %f %f %.1f\n",i/10.0+0.05,statang[i]/countnum[i],statcosang[i]/countnum[i],statsz[i]/countnum[i],countnum[i]);
	}
	fprintf(fp,"\n---for the interfacial waters---\n\n");
	fprintf(fp,"z avg-ang avg-cosang avg-sz\n\n");
	for(i=0;i<=Nbins;i++)
	{
		if(icount[i] == 0)
			fprintf(fp,"%.3f %f %f %f %.1f\n",i/10.0+0.05,0.0,0.0,0.0,0.0);   // change for diff. z resolution
		else if(icount[i] > 0)
			fprintf(fp,"%.3f %f %f %f %.1f\n",i/10.0+0.05,iang[i]/countnum[i],icosang[i]/icount[i],isz[i]/icount[i],icount[i]);
	}
	fprintf(fp,"\n---for the bulk-like waters---\n\n");
	fprintf(fp,"z avg-ang avg-cosang avg-sz\n\n");
	for(i=0;i<=Nbins;i++)
	{
		if(bcount[i] == 0)
			fprintf(fp,"%.3f %f %f %f %.1f\n",i/10.0+0.05,0.0,0.0,0.0,0.0);   // change for diff. z resolution
		else if(bcount[i] > 0)
			fprintf(fp,"%.3f %f %f %f %.1f\n",i/10.0+0.05,bang[i]/bcount[i],bcosang[i]/bcount[i],bsz[i]/bcount[i],bcount[i]);
	}



	fclose(fip);
	fclose(fp);
	return(0);
}


