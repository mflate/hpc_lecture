/******************************************************************************
 * Copyright (c) 2017, NVIDIA CORPORATION.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

#pragma once

/**
 * \file
 * Epilogue operation to compute final output
 */

namespace cutlass {
namespace gemm {

    //// Used by GEMM to compute the final result C <= alpha * accumulator + beta * C
    class blas_scaled_epilogue
    {
    public:

        float alpha;
        float beta;

        inline __device__ __host__
        blas_scaled_epilogue(
            float alpha,
            float beta)
        :
            alpha(alpha),
            beta(beta)
        {}


        /// Epilogue operator
        inline __device__ __host__
        float operator()(
            float accumulator,
            float c,
            size_t idx) const
        {
            return float(alpha * float(accumulator) + beta * float(c));
        }

        /**
         * Configure epilogue as to whether the thread block is a secondary
         * accumulator in an inter-block k-splitting scheme
         */
        inline __device__
        void set_secondary_accumulator()
        {
            beta = float(1);
        }


        /// Return whether the beta-scaled addend needs initialization
        inline __device__
        bool must_init_addend()
        {
            return (beta != float(0));
        }
    };




} // namespace gemm
} // namespace cutlass
