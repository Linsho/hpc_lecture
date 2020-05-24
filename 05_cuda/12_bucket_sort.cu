#include <cstdio>
#include <cstdlib>
#include <vector>


__global__ void init(int *a, int n){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= n) return;
  a[i] = 0;
}

__global__ void packing(int *key, int *bucket, int n){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= n) return;
  atomicAdd(&bucket[key[i]], 1);
} 

__global__ void scan(int *a, int *b, int n){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= n) return;
  for(int j = 1;j < n;j <<= 1){
    b[i] = a[i];
    __syncthreads();
    if(i-j >= 0) a[i] += b[i-j];
    __syncthreads();
  }
}

__global__ void unpacking(int *key, int *bucket, int n, int range){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= n) return;
  int top = range-1;
  int bottom = -1;
  int middle;
  while(top-bottom > 1){
    middle = (top+bottom)/2;
    if(i >= bucket[middle]) bottom = middle;
    else top = middle;
  }
  key[i] = top;
}  

int main() {
  int n = 50;
  const int m = 1024;
  int range = 5;
  int *key;
  cudaMallocManaged(&key, n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
  int *bucket;
  cudaMallocManaged(&bucket, range*sizeof(int));
  init<<<(n+m-1)/m, m>>>(bucket, n);
  packing<<<(n+m-1)/m, m>>>(key, bucket, n);
  int *scan_mem;
  cudaMallocManaged(&scan_mem, range*sizeof(int));
  scan<<<(range+m-1)/m, m>>>(bucket, scan_mem, range);
  unpacking<<<(n+m-1)/m, m>>>(key, bucket, n, range);  
  cudaDeviceSynchronize();
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
