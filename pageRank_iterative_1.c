#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scan_q1.h"
//#include <omp.h>

int n_vertex;
void read_from_file(int num_t, int *row_ar, int* col_ar ,FILE* f,int nelem)
{

	while(!feof(f))
	{
		
		int row,col;
		fscanf(f,"%d",&row);
		char c=getc(f);
		while((c=getc(f))!='\n'&& !(feof(f)))
		{
			fscanf(f,"%d",&col);
			col_ar[col]+=1;
			//m->elements[row-1][col-1]=1.0;
		}
			
	}
	int argv[2]; /*Pass as an argument in main function*/
	argv[0]=n_vertex+1;
	argv[1]=num_t;
	//printf("Columns\n");
	int i=0;
		
		
	scan(argv,col_ar);
	int tmp_col_ar[n_vertex+1];

	for(i=0;i<=n_vertex;i++)
		tmp_col_ar[i]=col_ar[i];	

	/*Storing the row_values of each element column-wise*/
	
	/*printf("Columns\n");
	for(i=0;i<=n_vertex;i++)
		printf("%d\n",col_ar[i]);
	printf("Columns\n");*/
		FILE* f1= fopen("testfile_links.txt","r");

	while(!feof(f1))
	{
		int row,col;
		fscanf(f1,"%d",&row);
		char c;
		while((c=getc(f1))!='\n'&& !(feof(f1)))
		{
			fscanf(f1,"%d",&col);
			row_ar[tmp_col_ar[col-1]]=row-1; /*Rows numbering start from zero*/
			tmp_col_ar[col-1]+=1;
			//m->elements[row-1][col-1]=1.0;
		}
			
	}
	/*printf("Row aray\n");
	for(i=0;i<nelem;i++)
		printf("%d\n",row_ar[i]);*/
	//f1= fopen("testfile_links.txt","r");
	fclose(f1);

}

	
int wordcount_from_file(int* out_deg,FILE* f)
{
	int count=0,i=0;
	while(!feof(f))
	{
		
		int row,col;
		fscanf(f,"%d",&row);
		char c=getc(f);
		int row_count=0;
		while((c=getc(f))!='\n'&& !(feof(f)))
		{
			fscanf(f,"%d",&col);
			++row_count;
			//m->elements[row-1][col-1]=1.0;
		}
		count+=row_count;
		if(row_count==0)
			row_count = n_vertex;
		if(i==n_vertex)break;

		out_deg[i]=row_count;
		i++;
			
	}

	//printf("%d\n",i);
	return count;
}


void initialize_vector(float* vector)
{
	int i=0;
	for(i=0;i<n_vertex;i++)	
	{
		//vector[i]=1.0/n_vertex;
		//printf("%lf",vector[i]);
		vector[i]=1.0/n_vertex;
	}
	
}

void vector_normalize(float* vector)
{
    int i=0;
    float sum = 0.0;
//    #pragma omp parallel for reduction(+:sum)
    for (i = 0; i < n_vertex; ++i)
    {
	sum += vector[i];
    }
     
//     #pragma omp parallel for 
    for (i = 0; i < n_vertex; ++i)
    {
	vector[i] /= sum;
    }
}


int main(int argc, char** argv)
{
	int num_t=atoi(argv[2]);
	n_vertex=atoi(argv[1]);
	int out_deg[n_vertex]; //Store the out-degree of each vertex
	float vector[n_vertex];
	float tmp_vector[n_vertex]; 
	
	FILE* f= fopen("testfile_links.txt","r");
	int row_ar_size=wordcount_from_file(out_deg,f);
	
	int i=0;	
	/*for(i=0;i<n_vertex;i++)
	{
		printf("%d\n",out_deg[i]);		
	}*/

	printf("Num elem %d\n",row_ar_size);
	int col_val[n_vertex+1];
	int row_val[row_ar_size];
	for(i=0;i<=n_vertex;i++)
		col_val[i]=0;	
	f= fopen("testfile_links.txt","r");
	
	read_from_file(num_t, row_val, col_val ,f,row_ar_size); //store the matrix in compressed row format
	
	int j=0,k=0;
	int niter=14;
	initialize_vector(vector);
	float sum=0.0,imp=0.0;
	float p=1;
	for(j=0;j<niter;j++)
	{
		sum=0.0; 
	//	#pragma omp parallel num_threads(num_t) 
	//	{	
	//		#pragma omp for reduction(+:sum)
			for(k=0;k<n_vertex;k++)
				sum+=vector[k];
	//		#pragma omp barrier
	//		#pragma omp parallel for private(imp)
			for(k=0;k<n_vertex;k++)
			{
				int start = col_val[k], end=col_val[k+1],m;
				imp=0.0;
				for(m=start;m<end;m++)
				{
					imp=imp+vector[row_val[m]]/out_deg[row_val[m]];
				}
				tmp_vector[k]=imp*p+((1-p)*(sum/n_vertex));
			
			}
	//		#pragma omp barrier
			
	//		#pragma omp parallel for
			for(k=0;k<n_vertex;k++)
			{
				vector[k]=tmp_vector[k]; /*copying the data back to the original vector*/
			}

			
	//	}
		vector_normalize(vector);
	}
	printf("PageRank\n");
	
sum=0;
for(k=0;k<n_vertex;k++)
{
		printf("%lf\n", vector[k]);
		sum+=vector[k];

}
printf("sum\n");
printf("%lf\n",sum);
	return 0;
}
