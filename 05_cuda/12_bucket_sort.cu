#include <cstdio>
#include <cstdlib>
#include <vector>

#include <cooperative_groups.h>
using namespace cooperative_groups;

// To run:
// nvcc 12_bucket_sort.cu -arch=sm_60 -rdc=true -L/usr/lib/x86_64-linux-gnu 

__global__ void initialize(int* bucket) {
  bucket[threadIdx.x] = 0;
}

__device__ int warpSum(int a) {
for (int offset=16; offset>0; offset >>= 1)
   a += __shfl_down_sync(0xffffffff, a, offset);
   return a;
}

// A warpsum for every bucket, then every thread adds
__global__ void sort(int*bucket, int *key, int*keymask, int range, int n) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  grid_group grid = this_grid();
  if (i >= n) return;
  keymask[i] = 0;
  for (int j = 0; j < range; j++) {
    if (key[i] == j)    keymask[i] = 1;
    grid.sync();
    int b = warpSum(keymask[i]); keymask[i] = 0;
    if ((i & 31) == 0) 
      atomicAdd(&bucket[j], b);    
    grid.sync();
  } 
  
  int sum = 0;
  for (int k = 0; k < range; k++) {
    sum+=bucket[k];
    if (i < sum) {key[i] = k; return;}
  }
}


int main() {
  int M = 32;
  int n = 50;
  int range = 5;
  // Changed datatype of key to make compatible with cuda
  int *key;
  cudaMallocManaged(&key,n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  int *bucket;
  cudaMallocManaged(&bucket, range*sizeof(int));
  initialize<<<1,range>>>(bucket);
  
  int *keymask;
  cudaMallocManaged(&keymask,n*sizeof(int));

  void *args[] = {(void *)&bucket, (void *)&key, (void *)&keymask, (void *)&range, (void *)&n};
  cudaLaunchCooperativeKernel((void*)sort,(n+M-1)/M, M, args);
  
  cudaDeviceSynchronize();
  printf("\n");
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
