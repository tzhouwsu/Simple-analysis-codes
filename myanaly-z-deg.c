#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

// here I update the code to get the density and degree histogram within each Zbins for interfacial waters and other waters
//#define Zbins 102  // z resolution of 1 Ang with 100 Ang in Z total
#define Zbins 1002  // Z resolution of 0.1 Ang within 100 Ang in Z, from -50 to 50
#define FLN 1000
#define Nnodes 1204  // the number of water molecules
#define Nsnaps 20000    // number of snapshots
#define Gibbs 33.06   // the Gibbs dividing surface, to analyze 2-Ang layering, from -6, -4, ... 8, 10, 12


int getdata(float ave[4], int snap, float numi[Zbins+1], float numb[Zbins+1], int degsi[Zbins+1][10], int degsb[Zbins+1][10])
{
	FILE *fp1,*fp2,*fp3;
	char filename[FLN];
	sprintf(filename,"interfacial-id-%d",snap);  // interfacial water ids
	if((fp1=fopen(filename,"r"))==NULL)
	{
		printf("Can not find interfacial-id file at snap %d\n",snap);
		return(-1);
	}
	sprintf(filename,"water%d.xyz",snap);  // xyz coordinate file
	if((fp2=fopen(filename,"r"))==NULL)
	{
		printf("Can not find water xyz file at snap %d\n",snap);
		fclose(fp1);
		return(-2);
	}
	sprintf(filename,"gg.%d",snap);  // water-water graphgeod files, corrected
	if((fp3=fopen(filename,"r"))==NULL)
	{
		printf("Can not find gg.%d file at snap %d\n",snap,snap);
		fclose(fp1);
		fclose(fp2);
		return(-3);
	}
	
	int i,j,id;
	char buffer[FLN];
	int index[Nnodes+1]; // labeling whether the water is interacting with hexane

	for(i=0;i<=Nnodes;i++)
		index[i]=0;

	rewind(fp1);
	while( fscanf(fp1,"%d",&id)==1 )
	{
//		printf(" %d %d\n",snap,id);
		index[id]=1;   // this is an interfacial water
		if(fgets(buffer,sizeof(buffer),fp1) == NULL)
			break;
	}
	fclose(fp1);


	float Coord[Nnodes][9];
	char lab;
	float x,y,z;
	int rderr,line; // read error
	for(i=0;i<=Nnodes;i++)
		for(j=0;j<9;j++)
			Coord[i][j]=-1.0;

	rewind(fp2);
	fgets(buffer,sizeof(buffer),fp2);
	fgets(buffer,sizeof(buffer),fp2);
	rderr=0;line=0;
	while(fscanf(fp2,"%c %e %e %e",&lab,&x,&y,&z)==4)
	{
		fgets(buffer,sizeof(buffer),fp2);
		line++;
		if(lab == 'O' && line>=1 && line<=Nnodes)
		{
			Coord[line][0]=x;   // the coordinates of oxygen
			Coord[line][1]=y;
			Coord[line][2]=z;

			fscanf(fp2,"%c %e %e %e",&lab,&x,&y,&z);
			fgets(buffer,sizeof(buffer),fp2);
			Coord[line][3]=x;   // the coordinates of hydrogens
			Coord[line][4]=y;
			Coord[line][5]=z;

			fscanf(fp2,"%c %e %e %e",&lab,&x,&y,&z);
			fgets(buffer,sizeof(buffer),fp2);
			Coord[line][6]=x;
			Coord[line][7]=y;
			Coord[line][8]=z;

		}
		else
		{
			rderr = 1;   // if there is error in reading the xyz file
			break;
		}
	}
	fclose(fp2);
	if(rderr==1)
	{
		fclose(fp3);
		printf("Error in reading the xyz file at snap %d\n",snap);
		return(-4);
	}


	int HBij[Nnodes+1][Nnodes+1];  // the water-water HBs at this snapshot
	int wati,watj,atmi,atmj;
	for(wati=0;wati<=Nnodes;wati++)
		for(watj=0;watj<=Nnodes;watj++)
			HBij[wati][watj]=0;

	rewind(fp3);
	while(fscanf(fp3,"%d %d %d %d",&wati,&watj,&atmi,&atmj)==4)
	{
		fgets(buffer,sizeof(buffer),fp3);  // skip the atom label
		if(wati>=1 && wati<=Nnodes && watj>=1 && watj<=Nnodes)  // make sure these are water molecules
		{
			if(atmi==1) // water i is the acceptor, contributing Oxygen
				HBij[watj][wati]=1;    // indicate there is a HB  watj --> wati
			if(atmj==1)
				HBij[wati][watj]=1;    // indicate there is a HB  wati --> watj
		}
	}
	fclose(fp3);


	float sumzl,sumnuml,sumzr,sumnumr;
	int bin,degcount;
	float shiftz;   // the shiftted z coordinate with respect to the Gibbs dividing surface

	sumzl=sumnuml=sumzr=sumnumr=0.0;
	for(i=1;i<=Nnodes;i++)
	{

		// this is for 1 Ang or 0.1 Ang layering analysis
		z = Coord[i][2];   // the Z coordinates of Oxygen of the water
		bin = (int) ((z+50)*10);     // 0.1 Ang resolution, z from -50 to 50 
//		bin = (int) (z+50);     // 1 Ang resolution, z from -50 to 50 
		if(bin<0)
			bin=0;
		else if(bin>Zbins)
			bin=Zbins;

/*
		// this is for 2 Ang layring analysis
		z = Coord[i][2];
		bin = 0;
		shiftz = z - Gibbs;  // shiftted z coordinate
		if(shiftz >= -7 && shiftz < -5)   // get statistics of 2 Ang layering results
			bin=1;
		if(shiftz >= -5 && shiftz < -3)
			bin=2;
		if(shiftz >= -3 && shiftz < -1)
			bin=3;
		if(shiftz >= -1 && shiftz < 1)
			bin=4;
		if(shiftz >= 1 && shiftz < 3)
			bin=5;
		if(shiftz >= 3 && shiftz < 5)
			bin=6;
		if(shiftz >= 5 && shiftz < 7)
			bin=7;
		if(shiftz >= 7 && shiftz < 9)
			bin=8;
		if(shiftz >= 9 && shiftz < 11)
			bin=9;
		if(shiftz >= 11 && shiftz < 13)
			bin=10;
*/
	
  	//	printf(" %d %d %f %f %f\n",snap,i,Coord[i][0],Coord[i][1],Coord[i][2]);

		if(index[i] == 1)    // this is the interfacial water
		{
			numi[bin] += 1;    // count the number of interfacial water in a Z-bin range

			if(z<0.0)     // the left water-part of box, the box is from -50 to 50
			{
				sumnuml += 1.0;   // the total number of interfacial water at the left interface
				sumzl += z;
			}
			else
			{
				sumnumr += 1.0;   // the total number of interfacial water at the right interface
				sumzr += z;
			}

			// here I want to calculate the degree histogram, at each bin layer
			degcount = 0; // initial the degree count of water i
			for(j=1;j<=Nnodes;j++)
			{
				if(HBij[i][j]==1)
					degcount += 1;
				if(HBij[j][i]==1)
					degcount += 1;
			}
			degsi[bin][degcount] += 1;   // degcount is the degree of water i

		}
		else //not interfacial water
		{
			numb[bin] += 1;  // count the number of non-interfacial water in a Z-bin range

			// here I want to calculate the degree histogram, at each bin layer
			degcount = 0; // initial the degree count of water i
			for(j=1;j<=Nnodes;j++)
			{
				if(HBij[i][j]==1)
					degcount += 1;
				if(HBij[j][i]==1)
					degcount += 1;
			}
			degsb[bin][degcount] += 1;   // degcount is the degree of water i
	
		}

	}

	ave[0] += sumzl;
	ave[1] += sumnuml;
	ave[2] += sumzr;
	ave[3] += sumnumr;

	return(0);
}

int main()
{
	FILE *fp,*fop;
	fp = fopen("result-density-z","w");   // the output file
	fop = fopen("result-degree-z","w");   // the degree histogram

	float average[4];    // to calculate the average z position
	float waternumi[Zbins+1],waternumb[Zbins+1];  // the number density of interfacial water and others, as a function of z
	int degint[Zbins+1][10],degbulk[Zbins+1][10];  // the degree histogram of interfacial water and others, as a function of z
	int t,ideg,sym,degsumi,totsumi,degsumb,totsumb;

	for(t=0;t<=Zbins;t++)
		waternumi[t]=waternumb[t]=0.0;
	for(t=0;t<=Zbins;t++)
		for(ideg=0;ideg<10;ideg++)
			degint[t][ideg]=degbulk[Zbins+1][ideg]=0;

	sym=0;
	for(t=1;t<=Nsnaps;t++)
	{
		printf("for snapshot %d\n",t);
		sym = getdata(average,t,waternumi,waternumb,degint,degbulk);
		if(sym!=0)
		{
			printf("Error in func 'getdata', at snap %d\n",t);
			break;
		}
	}


	fprintf(fp,"\nZ-coord count-int\n\n");
	for(t=0;t<=Zbins;t++)
		fprintf(fp,"%.2f %.2f\n",(t)/10.0-49.95,waternumi[t]);   // 0.1 Ang resolution
//		fprintf(fp,"%.1f %.1f\n",(t-50.0)/1.0+0.5,waternumi[t]);   // 1 Ang resolution
	fprintf(fp,"\n\n");
	fprintf(fp,"\nZ-coord count-bulk\n\n");
	for(t=0;t<=Zbins;t++)
		fprintf(fp,"%.2f %.2f\n",(t)/10.0-49.95,waternumb[t]);   // 0.1 Ang resolution
//		fprintf(fp,"%.1f %.1f\n",(t-50.0)/1.0+0.5,waternumb[t]);   // 1 Ang resolution
	fprintf(fp,"\n\n");

	if(average[1]>0 && average[3]>0)   // this is used to get the average Z position for the interfacial water
		fprintf(fp,"%f %f %f : %f %f %f\n",average[0],average[1],average[0]/average[1],average[2],average[3],average[2]/average[3]);
	else
		fprintf(fp,"%f %f : %f %f\n",average[0],average[1],average[2],average[3]);
	fprintf(fp,"\n\n");



	fprintf(fop,"\nZ-coord count-degrees:(0 to 9) for interfacial water\n");
	for(t=0;t<=Zbins;t++)   // print the degree histogram for interfacial water
	{
		degsumi=totsumi=0;
		for(ideg=0;ideg<10;ideg++)
		{
			degsumi += degint[t][ideg];
			totsumi += ideg * degint[t][ideg];
		}

		if(degsumi > 0)
			fprintf(fop,"%.2f  %f : %d %d %d %d %d %d %d %d %d %d\n",(t)/10.0-49.95,(totsumi+0.0)/(degsumi+0.0),degint[t][0],degint[t][1],degint[t][2],degint[t][3],degint[t][4],degint[t][5],degint[t][6],degint[t][7],degint[t][8],degint[t][9]);  // 0.1 Ang
//			fprintf(fop,"%.1f  %f : %d %d %d %d %d %d %d %d %d %d\n",(t-50.0)/1.0+0.5,(totsumi+0.0)/(degsumi+0.0),degint[t][0],degint[t][1],degint[t][2],degint[t][3],degint[t][4],degint[t][5],degint[t][6],degint[t][7],degint[t][8],degint[t][9]);  // 1 Ang
		else
			fprintf(fop,"%.2f %f : %d %d %d %d %d %d %d %d %d %d\n",(t)/10.0-49.95,0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);  // 0.1 Ang
//			fprintf(fop,"%.1f %f : %d %d %d %d %d %d %d %d %d %d\n",(t-50.0)/1.0+0.5,0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);  // 1 Ang
	}
	fprintf(fop,"\n\n");


	fprintf(fop,"\nZ-coord count-degrees:(0 to 9) for non-interfacial water\n");
	for(t=0;t<=Zbins;t++)   // print the degree histogram for interfacial water
	{
		degsumb=totsumb=0;
		for(ideg=0;ideg<10;ideg++)
		{
			degsumb += degbulk[t][ideg];
			totsumb += ideg * degbulk[t][ideg];
		}

		if(degsumb > 0)
			fprintf(fop,"%.2f  %f : %d %d %d %d %d %d %d %d %d %d\n",(t)/10.0-49.95,(totsumb+0.0)/(degsumb+0.0),degbulk[t][0],degbulk[t][1],degbulk[t][2],degbulk[t][3],degbulk[t][4],degbulk[t][5],degbulk[t][6],degbulk[t][7],degbulk[t][8],degbulk[t][9]);  // 0.1 Ang
//			fprintf(fop,"%.1f  %f : %d %d %d %d %d %d %d %d %d %d\n",(t-50.0)/1.0+0.5,(totsumb+0.0)/(degsumb+0.0),degbulk[t][0],degbulk[t][1],degbulk[t][2],degbulk[t][3],degbulk[t][4],degbulk[t][5],degbulk[t][6],degbulk[t][7],degbulk[t][8],degbulk[t][9]); // 1 Ang
		else
			fprintf(fop,"%.2f %f : %d %d %d %d %d %d %d %d %d %d\n",(t)/10.0-49.95,0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);  // 0.1 Ang
//			fprintf(fop,"%.1f %f : %d %d %d %d %d %d %d %d %d %d\n",(t-50.0)/1.0+0.5,0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);  // 1 Ang
	}
	fprintf(fop,"\n\n");



	fclose(fp);
	fclose(fop);
	return(0);

}




