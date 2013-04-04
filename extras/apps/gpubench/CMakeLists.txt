cmake_minimum_required (VERSION 2.8)
project (sandbox_weese_apps_gpubench)

#find_package (OpenCL)
#if (OPENCL_FOUND)
#  add_definitions(-DSEQAN_HAS_OPENCL=1)
#  include_directories(${OPENCL_INCLUDE_DIR})
#endif (OPENCL_FOUND)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../include)


# Update the list of file names below if you add source files to your test.
seqan_add_cuda_executable(gpubench      gpubench.cu gpubench.h myers.h testset.h)
#seqan_add_cuda_test_executable(gpubench      gpubench.cu myers.h testset.h)
#seqan_add_cuda_test_executable(gpubench_host gpubench.cu myers.h testset.h)

#if (OPENCL_FOUND)
		#	target_link_libraries(gpu_myers ${OPENCL_LIBRARIES})
#endif (OPENCL_FOUND)
