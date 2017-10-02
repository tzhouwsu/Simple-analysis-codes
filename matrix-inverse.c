
// here I want to test Gauss-Jordan method of calculating matrix inverse

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#define Ndummy 100

int myprint(float **array, int num);

int main(int argc, char *argv[])
{
	FILE *fip;

	if(argc != 2)
	{
		printf("Error, in using this code\n");
		exit(-1);
	}

	if( (fip = fopen(argv[1],"r")) == NULL )
	{
		printf("Error: cannot find file %s\n",argv[1]);
		exit(-1);
	}

	int row,column,irow,icolumn,m,jrow;
	char buffer[100];
	char *token;
	float temp;
	float dummyrow[Ndummy];
	float **GJ_matrix; // this is a GJ_matrix (initially, original matrix with identify matrix)

	rewind(fip);
	row = column = irow = icolumn = 0;
	while(fgets(buffer,sizeof(buffer),fip) != NULL)
	{
		token = strtok(buffer," \n");
		if(token == NULL)
			continue;
		else
		{
			row++;
			temp = atof(token);
			icolumn = 1;
			dummyrow[icolumn] = temp;
			while( (token = strtok(NULL," \n")) != NULL )
			{
				temp = atof(token);
				icolumn++;
				dummyrow[icolumn] = temp;
			}

			if(icolumn >= Ndummy)
			{
				printf("error, size of dummy array not enough, change 'Ndummy' in source code\n");
				break;
			}

			if(row == 1)
			{
				column = icolumn;
				GJ_matrix = (float **) malloc(column * sizeof(float *)); // assume square matrix, N rows
				for(m=0;m<column;m++)
					GJ_matrix[m] = (float *) calloc(column*2, sizeof(float)); // 2N columns

				for(m=1;m<=icolumn;m++)
					GJ_matrix[row-1][m-1] = dummyrow[m];
				GJ_matrix[row-1][icolumn-1+row] = 1.0;  // this is to put in the identify matrix
			}
			else
			{
				if(row > column)
				{
					printf("Error: row > column\n");
					break;
				}

				if(icolumn != column)
				{
					printf("Error: not reading a matrix, row %d\n",row);
					break;
				}
				else
				{
					for(m=1;m<=icolumn;m++)
						GJ_matrix[row-1][m-1] = dummyrow[m];
					GJ_matrix[row-1][icolumn-1+row] = 1.0;  // this is to put in the identify matrix
				}
			}
		}
	}
	fclose(fip);

	if(row != column)
	{
		printf("Error: this is not square matrix %d * %d\n",row,column);
		return(-1);
	}


	myprint(GJ_matrix,column);

	// below I'm using Gauss-Jordan method to calculate the inverse of the matrix
	for(irow=0; irow<row; irow++)
	{
		if(GJ_matrix[irow][irow] == 0)  // if this row_i start with 0, find another row_j, and swap them
		{
			for(jrow=irow+1; jrow<row; jrow++)
				if(GJ_matrix[jrow][irow] != 0)
					break;
			if(jrow >= row)
			{
				printf("error, matrix invertible\n");
				break;
			}
			else
			{
				for(icolumn =0; icolumn<column*2; icolumn++)
				{
					temp = GJ_matrix[irow][icolumn];
					GJ_matrix[irow][icolumn] = GJ_matrix[jrow][icolumn];
					GJ_matrix[jrow][icolumn] = temp;
				}
			}

		} // after this for loop, it is ensured that the first element of row_i is not zero

		temp = GJ_matrix[irow][irow];
		for(icolumn=0; icolumn<column*2; icolumn++)
			GJ_matrix[irow][icolumn] = GJ_matrix[irow][icolumn] / temp; // make element [row_i,column_i] to 1

//	myprint(GJ_matrix,column);

		for(jrow=0; jrow<row;jrow++) // this is to make all other rows zero in column_i, except for row_i
		{
			if(jrow == irow)
				continue;
			else
			{
				temp = GJ_matrix[jrow][irow];
				for(icolumn=0; icolumn<column*2; icolumn++)
					GJ_matrix[jrow][icolumn] = GJ_matrix[jrow][icolumn] - GJ_matrix[irow][icolumn] * temp; // null the jrow in column_i
			}
		}

//	myprint(GJ_matrix,column);

	}
	

	myprint(GJ_matrix,column);

	return(0);
}


int myprint(float **array, int num)
{
	int i,j;

	printf("\n");
	for(i=0;i<num;i++)
	{
		for(j=0;j<num*2;j++)
			printf("%f ",array[i][j]);
		printf("\n");
	}
	printf("\n");

	return(0);
}

