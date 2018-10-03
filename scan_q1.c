#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>


int* scan(int* argv, int* a)
{
  int n = argv[0];
  int num_t = argv[1];
  int logn = (int)ceil(log(num_t)/log(2));
  int * scanSum = (int *) malloc(sizeof(int)*num_t);
  int i, j, k, start, end, sub_size, two_i, two_i_p1;
  double start_time_1, end_time_1;
  
  start_time_1 = omp_get_wtime();

  #pragma omp parallel num_threads(num_t)
  {
    #pragma omp for private(i, j, start, end)
      for(i = 0; i < num_t; i++){
        start = (i)*(n/num_t);
        end = (i+1)*(n/num_t);
        if(i == num_t-1)
          end = n;    
        for(j=start; j<end-1; j++)
          a[j+1] = a[j] + a[j+1];
        scanSum[i] = a[end-1];
      }
    #pragma omp barrier

    /*Up sweep*/
    #pragma omp for private(i, j, two_i, two_i_p1)    
      for(i = 0; i < logn; i++){
        two_i = 1 << i;
        two_i_p1 = 1 << (i+1);
        for(j = two_i - 1; j < (num_t - two_i); j+=two_i_p1)
          scanSum[j + two_i] = scanSum[j] + scanSum[j + two_i];
      }
    #pragma omp barrier

    /*Down sweep*/
    #pragma omp for private(i, j, two_i, two_i_p1)  
      for(i = logn; i > 0; i--){
        two_i = 1 << (i-1);
        two_i_p1 = 1 << (i-2);
        for(j = two_i - 1; j < (num_t - two_i_p1); j+=two_i)
          scanSum[j + two_i_p1] = scanSum[j] + scanSum[j + two_i_p1];
      }
    #pragma omp barrier

     #pragma omp master
     {
        for(i = num_t; i > 0; i--)
          	scanSum[i] = scanSum[i-1];
      	scanSum[0]=0;
  	  }
    #pragma omp barrier

    #pragma omp for private(i, j, start, end)
      for(i = 0; i < num_t; i++){
        start = (i)*(n/num_t);
        end = (i+1)*(n/num_t);
        if(i == num_t-1)
          end = n;
        for(j = start; j < end; j++)
          a[j] = a[j] + scanSum[i];
      } 
    #pragma omp barrier
  }
  end_time_1 = omp_get_wtime(); 

//  printf("%lf\n", end_time_1 - start_time_1); 

  return 0;
}


