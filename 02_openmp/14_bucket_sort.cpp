#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include <algorithm>

template<class T>
void prefix_sum(std::vector<T>& a) {
  int n = a.size();
  std::vector<T> b  = std::vector<T>(n);
#pragma omp parallel
  for(int j = 1;j < n;j <<= 1){
#pragma omp for
    for(int i = 0;i < n;i++){
      b[i] = a[i];
    }
#pragma omp for
    for(int i = j;i < n;i++){
      a[i] += b[i-j];
    }
  }
}

int main() {
  int n = 50;
  int range = 50;
  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");
<<<<<<< HEAD
  sort(key.begin(), key.end());
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  std::vector<int> bucket(range); 
#pragma omp parallel for
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
#pragma omp parallel for
  for (int i=0; i<n; i++) {
#pragma omp atomic update
    bucket[key[i]]++;
  }
  for (int i=0; i<range; i++) {
    printf("%d ", bucket[i]);
  }
  printf("\n");
  prefix_sum(bucket);//O(log(range))
#pragma omp parallel for
  for(int j = 0;j < n;j++){
    int i = upper_bound(bucket.begin(), bucket.end(), j)-bucket.begin();
    key[j] = i;
    // O(log(range))
=======

  std::vector<int> bucket(range,0); 
#pragma omp parallel for
  for (int i=0; i<n; i++)
#pragma omp atomic update
    bucket[key[i]]++;
  std::vector<int> offset(range,0);
  for (int i=1; i<range; i++) 
    offset[i] = offset[i-1] + bucket[i-1];
#pragma omp parallel for
  for (int i=0; i<range; i++) {
    int j = offset[i];
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
>>>>>>> 9c52e5565b22e57d279d08fe3117ffbea63ea44f
  }

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
