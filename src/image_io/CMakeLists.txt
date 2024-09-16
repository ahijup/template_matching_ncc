cmake_minimum_required (VERSION 3.8)

set(PROJECT_NAME image_io_module)

file(GLOB BMP_PATHS
    "LibBmp/*.h"
    "LibBmp/*.cpp"
)
file(GLOB JPEG_FILES
    "LibJpeg/*.h"
    "LibJpeg/*.c"
    )
file(GLOB ZLIB_FILES
    "Zlib/*.h"
    "Zlib/*.c"
    )
file(GLOB TIFF_FILES
    "LibTiff/*.h"
    "LibTiff/*.c"
    )
file(GLOB PNG_FILES
    "LibPng/*.h"
    "LibPng/*.c"
    )
add_library(${PROJECT_NAME}  
      ${BMP_PATHS}
      ${JPEG_FILES}
      ${ZLIB_FILES}
      ${TIFF_FILES}
      ${PNG_FILES}
      "image_api.cpp"
      "image_api.h"
)

target_link_libraries(${PROJECT_NAME} Shlwapi)
