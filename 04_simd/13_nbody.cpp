#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

/* The code turned out quite long, but counting the lines, there should be a significant
reduction in cycles, I estimated it was around 1/4 of the original. */

int main() {
  const int N = 8;
  // Adding array ind, in order to create a mask
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;

  }

  // Initializing reusable vectors
  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);
  __m256 mvec = _mm256_load_ps(m);

  __m256 fxvec = _mm256_setzero_ps();
  __m256 fyvec = _mm256_setzero_ps();
  
  __m256 zerovec = _mm256_setzero_ps();

  for(int i=0; i<N; i++) {
<<<<<<< HEAD

    // General N-body calculations
    __m256 xivec = _mm256_set1_ps(x[i]);
    __m256 yivec = _mm256_set1_ps(y[i]);
    __m256 rxvec = _mm256_sub_ps(xivec, xvec);
    __m256 ryvec = _mm256_sub_ps(yivec, yvec);
    __m256 rx2vec = _mm256_mul_ps(rxvec, rxvec);
    __m256 ry2vec = _mm256_mul_ps(ryvec, ryvec);
    __m256 r2vec = _mm256_add_ps(rx2vec, ry2vec);
    __m256 rrsqvec =  _mm256_rsqrt_ps(r2vec);
    __m256 rrsq2vec =  _mm256_mul_ps(rrsqvec,rrsqvec);
    __m256 rrsq3vec =  _mm256_mul_ps(rrsqvec,rrsq2vec);

    // Creating a mask for removing the same body
    __m256 mask = _mm256_cmp_ps(xvec, xivec, _CMP_NEQ_OQ);

    rrsq3vec = _mm256_blendv_ps(zerovec, rrsq3vec, mask);

    // More gerenal N-body calculations 
    __m256 tmpxvec =  _mm256_mul_ps(rxvec,mvec);
    __m256 fxivec = _mm256_mul_ps(tmpxvec,rrsq3vec);
    fxivec = _mm256_sub_ps(zerovec, fxivec);

    __m256 tmpyvec = _mm256_mul_ps(ryvec,mvec);
    __m256 fyivec = _mm256_mul_ps(tmpyvec,rrsq3vec);
    fyivec = _mm256_sub_ps(zerovec, fyivec);

    // Vector reduction
    __m256 fxipvec = _mm256_permute2f128_ps(fxivec,fxivec,1);
    fxivec = _mm256_add_ps(fxivec,fxipvec);
    fxivec = _mm256_hadd_ps(fxivec,fxivec);
    fxivec = _mm256_hadd_ps(fxivec,fxivec);

    __m256 fyipvec = _mm256_permute2f128_ps(fyivec,fyivec,1);
    fyivec = _mm256_add_ps(fyivec,fyipvec);
    fyivec = _mm256_hadd_ps(fyivec,fyivec);
    fyivec = _mm256_hadd_ps(fyivec,fyivec);

    // Adding the ith element to the fx and fy vectors
    fxvec = _mm256_blendv_ps(fxivec, fxvec, mask);
    fyvec = _mm256_blendv_ps(fyivec, fyvec, mask);
=======
    for(int j=0; j<N; j++) {
      if(i != j) {
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];
        float r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }
    }
    printf("%d %g %g\n",i,fx[i],fy[i]);
>>>>>>> 0294917e33255a08c281b7392dbbafc078a00cb3
  }
  _mm256_store_ps(fx, fxvec);
  _mm256_store_ps(fy, fyvec);

  // Printing the forces for debugging
  for(int i=0; i<N; i++) 
    printf("%d %f %f\n",i,fx[i],fy[i]); 
}
