#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scan_q1.h"
//#include <omp.h>
int n_vertex;

void read_from_file(int num_t,int *row_ar, int* col_ar ,FILE* f1,int nelem)
{

	int argv[2]; /*Pass as an argument to the scan function*/
	argv[0]=n_vertex+1;
	argv[1]=num_t;
	//printf("Columns\n");
		
		
	scan(argv,col_ar);
	int* tmp_col_ar=(int*)malloc((n_vertex+1)*sizeof(int));
	int i=0;
	
	for(i=0;i<=n_vertex;i++)
		tmp_col_ar[i]=col_ar[i];	

		//FILE* f1= fopen("testfile_links.txt","r");
	while(!feof(f1))
	{
		int row,col;
		fscanf(f1,"%d",&row);
		char c;
		while((c=getc(f1))!='\n'&& !(feof(f1)))
		{
			fscanf(f1,"%d",&col);
			row_ar[tmp_col_ar[col]]=row; /*Rows numbering start from zero*/
			tmp_col_ar[col]+=1;
			//m->elements[row-1][col-1]=1.0;
		}
			
	}
	
	/*printf("Row aray\n");
	for(i=0;i<nelem;i++)
		printf("%d\n",row_ar[i]);*/
	//f1= fopen("testfile_links.txt","r");
	fclose(f1);

}

	
int wordcount_from_file(int* out_deg,int* col_ar,FILE* f)
{
	int count=0;
	printf("Read from file\n");
	
	while(!feof(f))
	{
		
		int row,col;
		fscanf(f,"%d",&row);
		char c;
		int row_count=0;
		while((c=getc(f))!='-1'&& !(feof(f)))
		{
			fscanf(f,"%d",&col);
			//printf("col:%d \n",col);
			col_ar[col+1]+=1;
		
			++row_count;
			//m->elements[row-1][col-1]=1.0;
		}
		count+=row_count;
		//if(i==n_vertex)break;
	
		out_deg[row]=row_count;
		printf("row:%d out_deg:%d\n",row,out_deg[row]);
		//i++;
			
	}
printf("Read from file\n");
	
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
	int* out_deg=(int *)malloc(n_vertex*(sizeof(int))); //Store the out-degree of each vertex
	float* vector=(float *)malloc(n_vertex*(sizeof(int)));
	float* tmp_vector=(float *)malloc(n_vertex*(sizeof(int))); 
	int* col_val=(int *)malloc((n_vertex+1)*(sizeof(int)));
	int i=0;	
		
	for(i=0;i<=n_vertex;i++)
		col_val[i]=0;	
	
	for(i=0;i<=n_vertex;i++)
		out_deg[i]=n_vertex;	
		
	printf("check\n");
	FILE* f= fopen("search.txt","r");
	int row_ar_size=wordcount_from_file(out_deg,col_val,f);
	
	
	printf("Num elem %d\n",row_ar_size);
	int* row_val=(int *)malloc(row_ar_size*(sizeof(int)));
	f= fopen("search.txt","r");
	
	read_from_file(num_t,row_val, col_val ,f,row_ar_size); //store the matrix in compressed row format
	/*for(i=0;i<10;i++)
	{
		printf("%d\n",out_deg[i]);
		printf("col_val:%d\n",col_val[i+1]);		
	}*/
	
	int j=0,k=0;
	int niter=100;
	initialize_vector(vector);
	float sum=0.0,imp=0.0;
	float p=1;
	for(j=0;j<niter;j++)
	{
		sum=0.0; 
//		#pragma omp parallel num_threads(num_t) 
//		{	
//			#pragma omp for reduction(+:sum)
			for(k=0;k<n_vertex;k++)
				sum+=vector[k];
//			#pragma omp barrier
//			#pragma omp for private(imp)
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
//			#pragma omp barrier
			
//			#pragma omp for
			for(k=0;k<n_vertex;k++)
			{
				vector[k]=tmp_vector[k]; /*copying the data back to the original vector*/
			}

			
//		}
		vector_normalize(vector);
	}
	printf("PageRank\n");
	
sum=0;
/*
for(k=0;k<n_vertex;k++)
{
		printf("%lf\n", vector[k]);
		sum+=vector[k];

}*/
printf("sum\n");
printf("%lf\n",sum);

	return 0;
}
