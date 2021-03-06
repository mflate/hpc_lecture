#pragma once
#include <stdint.h>
#include "../util/util.h"
#include "block_task.h"
#include "k_split_control.h"
#include "epilogue_function.h"

namespace cutlass {
namespace gemm {


  __global__ void kernel(
                       int m,                      ///< Height in rows of op(A) and C
                       int n,                      ///< Width in columns of op(B) and C
                       int k,                      ///< Width in columns of op(A) and height in rows of op(B)
                       k_split_control k_split,    ///< Abstraction for controlling inter-block k-splitting
                       blas_scaled_epilogue op,           ///< Epilogue operation to update matrix C
                       float *d_a,               ///< Pointer to matrix A array values
                       float *d_b,               ///< Pointer to matrix B array values
                       float *d_c)               ///< Pointer to matrix C array values
{

    // Declare statically-allocated shared storage
    __shared__ typename block_task::scratch_storage_t smem;

    // Construct and run the task
    block_task(
        &smem,
        d_a,
        d_b,
        d_c,
        op,
        m,
        n,
        k,
        k_split).run();
}

/******************************************************************************
 * Dispatch stub
 ******************************************************************************/

/**
 * GEMM dispatch stub
 *
 * This function also serves as the autotuning entrypoint to evaluate different
 * tuning parameterizations of kernel.
 */
void dispatch(
    int             m,                              ///< Height in rows of op(A) and C
    int             n,                              ///< Width in columns of op(B) and C
    int             k,                              ///< Width in columns of op(A) and height in rows of op(B)
    float           alpha,
    float           beta,
    float         *d_a,                           ///< Device pointer to matrix A array values
    float         *d_b,                           ///< Device pointer to matrix B array values
    float         *d_c,                           ///< Device pointer to matrix C array values
    cudaStream_t    stream = 0)                     ///< CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.)      
                                                    ///  to check for errors.  Also causes launch configurations to be printed
                                                    ///  to the console if DEBUG is defined.  Default is \p false.
{
  
  blas_scaled_epilogue epilogue(alpha, beta);

  int BlockItemsX = 64;
  int BlockItemsY = 64;
  dim3 block = dim3(64);
  dim3 grid = dim3(
            (m + BlockItemsY - 1) / BlockItemsY,
            (n + BlockItemsX - 1) / BlockItemsX);
  int dynamic_smem_bytes = 0;
  int max_sm_occupancy = 8;
  int sm_count;
  get_sm_count(sm_count);
  int *d_flags;
  cudaGetSymbolAddress((void**) &d_flags, d_flags_split_k);

  k_split_control k_split(
                          d_flags,
                          sm_count,
                          max_sm_occupancy,
                          k,
                          8,
                          block,
                          grid);
  gemm::kernel
    <<< grid,
    block,
    dynamic_smem_bytes,
    stream >>>(
               m,
               n,
               k,
               k_split,
               epilogue,
               d_a,
               d_b,
               d_c);
}


} // namespace gemm
} // namespace cutlass
