
// this is to select the interfacial water molecules and the interfacial hexane molecules

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define Nwater 3205   // the number of water molecules
#define Nhexane 586   // the number of hexane molecules
#define Snapid 1000   // the snapshot id, this can be a user-specified number

int main()
{
	int i,j,k;
	float waterxyz[Nwater+1][9],hexanexyz[Nhexane+1][60]; // record the coordinates of waters and hexanes
	FILE *fgg,*fw,*fh,*fow1,*fow2,*foh1,*foh2;
	char filename[1000],buffer[1000];
	int id1,id2,nitfw,nitfh;
	int itfw[Nwater+1],itfh[Nhexane+1]; // record the interfacial water and interfacial hexane
	char lab;
	float tempx,tempy,tempz;
	char Molw[3],Molh[20];    // the atomic order of water and hexane

	Molw[0]='O';Molw[1]='H';Molw[2]='H';
	for(k=0;k<6;k++)
		Molh[k]='C';
	for(k=6;k<20;k++)
		Molh[k]='H';

	for(i=0;i<=Nwater;i++)
		itfw[i]=0;
	for(j=0;j<=Nhexane;j++)
		itfh[j]=0;

	for(i=0;i<=Nwater;i++)
		for(j=0;j<9;j++)
			waterxyz[i][j]=0.0;
	for(i=0;i<=Nhexane;i++)
		for(j=0;j<60;j++)
			hexanexyz[i][j]=0.0;
	
	// read the graphgeod file of water:hexane interactions, which are used to define "interfacial" molecules
	sprintf(filename,"Input.water%d.xyz.solB%d.xyz.GraphGeod",Snapid,Snapid);
	if((fgg=fopen(filename,"r"))==NULL)
	{
		printf("Error: cannot find the file 'Input.water%d.xyz.solB%d.xyz.GraphGeod'\n\n",Snapid,Snapid);
		exit(-1);
	}
	rewind(fgg);
	while(fscanf(fgg,"%d %d",&id1,&id2)==2)  // read the first two integer number, which are water-id and hexane-id
	{
		fgets(buffer,sizeof(buffer),fgg); // skip the rest of the line
		if( (id1>=1 && id1<=Nwater) && (id2>=Nwater+1 && id2<=Nwater+Nhexane) )
		{
			itfw[id1] = 1; // id1 is a water molecule
			itfh[id2-Nwater] = 2; // id2 is a hexane molecule
		}
		if( (id2>=1 && id2<=Nwater) && (id1>=Nwater+1 && id1<=Nwater+Nhexane) )
		{
			itfw[id2] = 1; // id2 is a water molecule
			itfh[id1-Nwater] = 2; // id1 is a hexane molecule
		}
	}
	fclose(fgg);

	// read the coordinates of water molecules
	sprintf(filename,"water%d.xyz",Snapid);
	if((fw=fopen(filename,"r"))==NULL)
	{
		printf("Error: cannot find the file 'water%d.xyz'\n\n",Snapid);
		exit(-2);
	}
	rewind(fw);
	fgets(buffer,sizeof(buffer),fw);
	fgets(buffer,sizeof(buffer),fw);
	i=0;
	while(fscanf(fw,"%c %e %e %e",&lab,&tempx,&tempy,&tempz)==4  && i<=Nwater)
	{
		fgets(buffer,sizeof(buffer),fw);
		if(lab=='O')
		{
			i++;
			k=1;
			waterxyz[i][0]=tempx;
			waterxyz[i][1]=tempy;
			waterxyz[i][2]=tempz;
			
			fscanf(fw,"%c %e %e %e",&lab,&tempx,&tempy,&tempz);   // read the first hydrogen atom
			fgets(buffer,sizeof(buffer),fw);
			if(lab=='H')
			{
				k++;
				waterxyz[i][3]=tempx;
				waterxyz[i][4]=tempy;
				waterxyz[i][5]=tempz;
			}

			fscanf(fw,"%c %e %e %e",&lab,&tempx,&tempy,&tempz);   // read the second hydrogen atom
			fgets(buffer,sizeof(buffer),fw);
			if(lab=='H')
			{
				k++;
				waterxyz[i][6]=tempx;
				waterxyz[i][7]=tempy;
				waterxyz[i][8]=tempz;
			}

//	printf(" %d %f %f %f %f %f %f %f %f %f\n",i,waterxyz[i][0],waterxyz[i][1],waterxyz[i][2],waterxyz[i][3],waterxyz[i][4],waterxyz[i][5],waterxyz[i][6],waterxyz[i][7],waterxyz[i][8]);
			if(k != 3)      // if the coodinate of oxygen is not followed with two hydrogens, skip the loop
				break;
		}
	}
	if(i != Nwater)
		printf("Warning: error in reading the file 'water%d.xyz'\n\n",Snapid);
	fclose(fw);

	// read the coodinates of hexane molecules
	sprintf(filename,"solB%d.xyz",Snapid);
	if((fh=fopen(filename,"r"))==NULL)
	{
		printf("Error: cannot find the file 'solB%d.xyz'\n\n",Snapid);
		exit(-3);
	}
	rewind(fh);
	fgets(buffer,sizeof(buffer),fh);
	fgets(buffer,sizeof(buffer),fh);
	j=0;
	while(fscanf(fh,"%c %e %e %e",&lab,&tempx,&tempy,&tempz)==4 && j<=Nhexane)
	{
		fgets(buffer,sizeof(buffer),fh);
		if(lab=='C')
		{
			j++;
			hexanexyz[j][0]=tempx;
			hexanexyz[j][1]=tempy;
			hexanexyz[j][2]=tempz;
//		printf(" %d %f %f %f\n",j,hexanexyz[j][0],hexanexyz[j][1],hexanexyz[j][2]);

			for(k=1;k<=19;)
			{
				fscanf(fh,"%c %e %e %e",&lab,&tempx,&tempy,&tempz);
				fgets(buffer,sizeof(buffer),fh);
				if(lab == Molh[k])
				{
					hexanexyz[j][k*3]=tempx;
					hexanexyz[j][k*3+1]=tempy;
					hexanexyz[j][k*3+2]=tempz;
					k=k+1;
				}
				else
					break;

			}
			if(k != 20)
				break;


		}

	}
	if(j != Nhexane)
		printf("Warning: error in reading the file 'solB%d.xyz'\n\n",Snapid);
	fclose(fh);

	// start to create the xyz files with interfacial water, interfacial hexane, other water, other hexane
	nitfw=0;
	for(i=1;i<=Nwater;i++)
		if(itfw[i] == 1)
			nitfw++;
	nitfh=0;
	for(j=1;j<=Nhexane;j++)
		if(itfh[j] == 2)
			nitfh++;

	sprintf(filename,"iwater%d.xyz",Snapid);
	fow1=fopen(filename,"w");
	sprintf(filename,"owater%d.xyz",Snapid);
	fow2=fopen(filename,"w");
	sprintf(filename,"ihexane%d.xyz",Snapid);
	foh1=fopen(filename,"w");
	sprintf(filename,"ohexane%d.xyz",Snapid);
	foh2=fopen(filename,"w");
	
	fprintf(fow1,"    %d\n\n",nitfw*3);
	fprintf(fow2,"    %d\n\n",(Nwater-nitfw)*3);
	fprintf(foh1,"    %d\n\n",nitfh*20);
	fprintf(foh2,"    %d\n\n",(Nhexane-nitfh)*20);
	for(i=1;i<=Nwater;i++)
	{
		if(itfw[i] == 1)
		{
			for(k=0;k<3;k++)
				fprintf(fow1," %c %f %f %f\n",Molw[k],waterxyz[i][k*3],waterxyz[i][k*3+1],waterxyz[i][k*3+2]);
		}
		else if(itfw[i] == 0)
		{
			for(k=0;k<3;k++)
				fprintf(fow2," %c %f %f %f\n",Molw[k],waterxyz[i][k*3],waterxyz[i][k*3+1],waterxyz[i][k*3+2]);
		}
	}
	for(j=1;j<=Nhexane;j++)
	{
		if(itfh[j] == 2)
		{
			for(k=0;k<20;k++)
				fprintf(foh1," %c %f %f %f\n",Molh[k],hexanexyz[j][k*3],hexanexyz[j][k*3+1],hexanexyz[j][k*3+2]);
		}
		else if(itfh[j] == 0)
		{
			for(k=0;k<20;k++)
				fprintf(foh2," %c %f %f %f\n",Molh[k],hexanexyz[j][k*3],hexanexyz[j][k*3+1],hexanexyz[j][k*3+2]);
		}
	}
	
	fclose(fow1);
	fclose(fow2);
	fclose(foh1);
	fclose(foh2);


	return(0);
}




