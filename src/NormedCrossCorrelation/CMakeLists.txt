find_package(OpenMP REQUIRED)

add_executable (NormedCrossCorrelation 
    "main.c" 
    "ImageData.h"
    "ncc.h"
    "ncc.c"
    "ncc_sum.h"
    "ncc_sum.c"
    "ncc_fft.h"
    "ncc_fft.c"
    "ncc_approx.h"
    "ncc_approx.c"
    "ncc_simd.h"
    "ncc_simd.c"
    "ncc_common.h"
    "ncc_common.c"
)

add_definitions(-DTEST_IMAGE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(NormedCrossCorrelation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../image_io)

target_link_libraries(NormedCrossCorrelation image_io_module OpenMP::OpenMP_CXX)
