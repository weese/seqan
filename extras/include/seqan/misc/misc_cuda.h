/*
 * Copyright (c) 2010-2011, NVIDIA Corporation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of NVIDIA Corporation nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef SEQAN_EXTRAS_MISC_CUDA_MISC_H_
#define SEQAN_EXTRAS_MISC_CUDA_MISC_H_

//#include <cuda_runtime.h>
//#include <thrust/version.h>
//#if THRUST_VERSION < 100600
//#include <thrust/detail/backend/cuda/arch.h>
//#else
//#include <thrust/system/cuda/detail/arch.h>
//#endif

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

/*
struct Arch
{
    static const __uint32 LOG_WARP_SIZE = 5;
    static const __uint32 WARP_SIZE     = 1u << LOG_WARP_SIZE;
};
*/

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function cudaPrintFreeMemory()
// --------------------------------------------------------------------------

inline void cudaPrintFreeMemory()
{
    size_t free, total;
    cudaMemGetInfo(&free, &total);

    std::cout << "Free" << free/1024/1024 <<  " of total " << total/1024/1024 << " MB\n";
}

// --------------------------------------------------------------------------
// Function cudaOccupancy()
// --------------------------------------------------------------------------

inline float cudaOccupancy()
{
#if defined(__CUDA_ARCH__)
    return __popc(__ballot(true)) / 32.0f;
#else
    return 1;
#endif
}

// --------------------------------------------------------------------------
// Function checkCudaError()
// --------------------------------------------------------------------------

//inline void checkCudaError(const char *message)
//{
//    cudaError_t error = cudaGetLastError();
//    if(error!=cudaSuccess) {
//        fprintf(stderr,"%s: %s\n", message, cudaGetErrorString(error) );
//        exit(1);
//    }
//}

/*
// granularity of shared memory allocation
inline size_t smem_allocation_unit(const cudaDeviceProp& properties)
{
    return 512;
}

// granularity of register allocation
inline size_t reg_allocation_unit(const cudaDeviceProp& properties)
{
    //return (properties.major < 2 && properties.minor < 2) ? 256 : 512;
    return
        (properties.major <= 1) ?
            (properties.minor <= 1 ? 256 : 512) :
            64;
}

// granularity of warp allocation
inline size_t warp_allocation_multiple(const cudaDeviceProp& properties)
{
    return 2;
}

inline size_t max_blocks_per_multiprocessor(const cudaDeviceProp& properties)
{
    return properties.major <= 2 ? 8 : 16;
}

inline size_t max_active_blocks_per_multiprocessor(const cudaDeviceProp&        properties,
                                                   const cudaFuncAttributes&    attributes,
                                                   size_t CTA_SIZE,
                                                   size_t dynamic_smem_bytes)
{
  // Determine the maximum number of CTAs that can be run simultaneously per SM
  // This is equivalent to the calculation done in the CUDA Occupancy Calculator spreadsheet
  const size_t regAllocationUnit      = reg_allocation_unit(properties);
  const size_t warpAllocationMultiple = warp_allocation_multiple(properties);
  const size_t smemAllocationUnit     = smem_allocation_unit(properties);
  const size_t maxThreadsPerSM        = properties.maxThreadsPerMultiProcessor;  // 768, 1024, 1536, etc.
  const size_t maxBlocksPerSM         = max_blocks_per_multiprocessor(properties);

  // Number of warps (round up to nearest whole multiple of warp size & warp allocation multiple)
  const size_t numWarps = util::round_i(util::divide_ri(CTA_SIZE, properties.warpSize), warpAllocationMultiple);

  // Number of regs is regs per thread times number of warps times warp size
  const size_t regsPerCTA = properties.major < 2 ?
      util::round_i(attributes.numRegs * properties.warpSize * numWarps, regAllocationUnit) :
      util::round_i(attributes.numRegs * properties.warpSize, regAllocationUnit) * numWarps;

  const size_t smemBytes  = attributes.sharedSizeBytes + dynamic_smem_bytes;
  const size_t smemPerCTA = util::round_i(smemBytes, smemAllocationUnit);

  const size_t ctaLimitRegs    = regsPerCTA > 0 ? properties.regsPerBlock      / regsPerCTA : maxBlocksPerSM;
  const size_t ctaLimitSMem    = smemPerCTA > 0 ? properties.sharedMemPerBlock / smemPerCTA : maxBlocksPerSM;
  const size_t ctaLimitThreads =                  maxThreadsPerSM              / CTA_SIZE;

  return _min( (__uint32)ctaLimitRegs, _min( (__uint32)ctaLimitSMem, _min((__uint32)ctaLimitThreads, (__uint32)maxBlocksPerSM)));
}

template <typename KernelFunction>
size_t max_active_blocks(KernelFunction kernel, const size_t CTA_SIZE, const size_t dynamic_smem_bytes)
{
    int            device;
    cudaGetDevice( &device );

    cudaDeviceProp properties;
    cudaGetDeviceProperties( &properties, device );

#ifdef __CUDACC__
    typedef void (*fun_ptr_type)();

    fun_ptr_type fun_ptr = reinterpret_cast<fun_ptr_type>(kernel);

    cudaFuncAttributes attributes;
    cudaFuncGetAttributes(&attributes, fun_ptr);

    return properties.multiProcessorCount * max_active_blocks_per_multiprocessor(properties, attributes, CTA_SIZE, dynamic_smem_bytes);
#else
    return 0u;
#endif
}

template <typename KernelFunction>
size_t num_registers(KernelFunction kernel)
{
#ifdef __CUDACC__
    typedef void (*fun_ptr_type)();

    fun_ptr_type fun_ptr = reinterpret_cast<fun_ptr_type>(kernel);

    cudaFuncAttributes attributes;
    cudaFuncGetAttributes(&attributes, fun_ptr);

    return attributes.numRegs;
#else
    return 0u;
#endif
}

template <typename KernelFunction>
size_t auto_blocksize(KernelFunction kernel, size_t dynamic_smem_bytes_per_thread = 0)
{
#if THRUST_VERSION < 100600
    return thrust::detail::backend::cuda::arch::max_blocksize_with_highest_occupancy(kernel, dynamic_smem_bytes_per_thread);
#else
    return thrust::system::cuda::detail::arch::max_blocksize_with_highest_occupancy(kernel, dynamic_smem_bytes_per_thread);
#endif
}

inline bool is_tcc_enabled()
{
    int            device;
    cudaDeviceProp device_properties;
    cudaGetDevice(&device);
    cudaGetDeviceProperties( &device_properties, device );
    return device_properties.tccDriver ? true : false;
}

struct cuda_error {};

inline void check_error(const char *message)
{
	cudaError_t error = cudaGetLastError();
	if(error!=cudaSuccess) {
		log_error(stderr,"%s: %s\n", message, cudaGetErrorString(error) );
        throw cuda_error();
	}
}

/// a generic syncthreads() implementation to synchronize contiguous
/// blocks of N threads at a time
///
template <__uint32 N>
SEQAN_FUNC
void syncThreads()
{
    #ifdef __CUDA_ARCH__
    if ((N > cuda::Arch::WARP_SIZE) || (is_pow2<N>() == false))
        __syncthreads();
    #endif
}
*/

}

#endif  // SEQAN_EXTRAS_MISC_CUDA_MISC_H_
