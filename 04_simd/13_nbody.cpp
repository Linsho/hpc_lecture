#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>
<<<<<<< HEAD

float _mm256_reduction(__m256 avec){
  __m256 bvec = _mm256_permute2f128_ps(avec,avec,1);
  bvec = _mm256_add_ps(bvec,avec);
  bvec = _mm256_hadd_ps(bvec,bvec);
  bvec = _mm256_hadd_ps(bvec,bvec);
  float a[8];
  _mm256_store_ps(a, bvec);
  return a[0]; 
}

void dbug_check(__m256 xvec){
  float x[8];
  _mm256_store_ps(x, xvec);
  for(int i = 0;i < 8;i++){
    printf("%g ", x[i]);
  }
  printf("\n");
}

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], idx[N];
=======

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
>>>>>>> 9c52e5565b22e57d279d08fe3117ffbea63ea44f
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
    idx[i] = (float)i;
  }
<<<<<<< HEAD
  for(int i=0; i<N; i++) {
    __m256 xvec = _mm256_load_ps(x);
    __m256 yvec = _mm256_load_ps(y);
    __m256 xivec = _mm256_set1_ps(x[i]);
    __m256 yivec = _mm256_set1_ps(y[i]);
    __m256 rxvec = _mm256_sub_ps(xivec, xvec);
    __m256 ryvec = _mm256_sub_ps(yivec, yvec);  
    __m256 rrvec = _mm256_add_ps(
      _mm256_mul_ps(rxvec, rxvec),
      _mm256_mul_ps(ryvec, ryvec)
    );
    __m256 one_rvec = _mm256_rsqrt_ps(rrvec);
    __m256 idxvec = _mm256_load_ps(idx);
    __m256 ivec = _mm256_set1_ps((float)i);
    __m256 mask = _mm256_cmp_ps(idxvec, ivec, _CMP_NEQ_US);
    __m256 zerovec = _mm256_set1_ps((float)0);
    __m256 mvec = _mm256_load_ps(m);
    __m256 fxvec = _mm256_mul_ps(rxvec, mvec); 
    fxvec = _mm256_mul_ps(fxvec, one_rvec);
    fxvec = _mm256_mul_ps(fxvec, one_rvec);
    fxvec = _mm256_mul_ps(fxvec, one_rvec);
    fxvec = _mm256_blendv_ps(zerovec, fxvec, mask);
    fx[i] -= _mm256_reduction(fxvec);
    
    __m256 fyvec = _mm256_mul_ps(ryvec, mvec); 
    fyvec = _mm256_mul_ps(fyvec, one_rvec);
    fyvec = _mm256_mul_ps(fyvec, one_rvec);
    fyvec = _mm256_mul_ps(fyvec, one_rvec);
    fyvec = _mm256_blendv_ps(zerovec, fyvec, mask);
    fy[i] -= _mm256_reduction(fyvec);
    
//    for(int j=0; j<N; j++) {
//      if(i != j) {
//        double rx = x[i] - x[j];
//        double ry = y[i] - y[j];
//        double r = std::sqrt(rx * rx + ry * ry);
//        fx[i] -= rx * m[j] / (r * r * r);
//        fy[i] -= ry * m[j] / (r * r * r);
//      }
//    }
    printf("%d %g %g\n",i,fx[i],fy[i]);
=======
  __m256 zero = _mm256_setzero_ps();
  for(int i=0; i<N; i+=8) {
    __m256 xi = _mm256_load_ps(x+i);
    __m256 yi = _mm256_load_ps(y+i);
    __m256 fxi = zero;
    __m256 fyi = zero;
    for(int j=0; j<N; j++) {
      __m256 dx = _mm256_set1_ps(x[j]);
      __m256 dy = _mm256_set1_ps(y[j]);
      __m256 mj = _mm256_set1_ps(m[j]);
      __m256 r2 = zero;
      dx = _mm256_sub_ps(xi, dx);
      dy = _mm256_sub_ps(yi, dy);
      r2 = _mm256_fmadd_ps(dx, dx, r2);
      r2 = _mm256_fmadd_ps(dy, dy, r2);
      __m256 mask = _mm256_cmp_ps(r2, zero, _CMP_GT_OQ);
      __m256 invR = _mm256_rsqrt_ps(r2);
      invR = _mm256_blendv_ps(zero, invR, mask);
      mj = _mm256_mul_ps(mj, invR);
      invR = _mm256_mul_ps(invR, invR);
      mj = _mm256_mul_ps(mj, invR);
      fxi = _mm256_fmadd_ps(dx, mj, fxi);
      fyi = _mm256_fmadd_ps(dy, mj, fyi);
    }
    _mm256_store_ps(fx+i, fxi);
    _mm256_store_ps(fy+i, fyi);
>>>>>>> 9c52e5565b22e57d279d08fe3117ffbea63ea44f
  }
  for(int i=0; i<N; i++)
    printf("%d %g %g\n",i,fx[i],fy[i]);
}
