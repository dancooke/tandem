cmake_minimum_required(VERSION 3.9)

include(CheckIPOSupported)

project(tandem)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(libdivsufsort)

add_library(libtandem STATIC tandem.hpp tandem.cpp)

set_target_properties(libtandem PROPERTIES OUTPUT_NAME "tandem")

set(WarningIgnores
    -Wno-unused-parameter
    -Wno-unused-function
    -Wno-missing-braces
    )

add_compile_options(-Wall -Wextra -Werror ${WarningIgnores})
target_link_libraries(libtandem libdivsufsort)
set_target_properties(libtandem PROPERTIES LINKER_LANGUAGE CXX)

#add_executable(tandem test.cpp ${TANDEM_SOURCES})
#target_link_libraries(tandem libtandem)
#install(TARGETS tandem DESTINATION ${tandem_SOURCE_DIR})

add_executable(tandem main.cpp ${TANDEM_SOURCES})
target_link_libraries(tandem libtandem)
install(TARGETS tandem DESTINATION ${tandem_SOURCE_DIR})

check_ipo_supported(RESULT ipo_supported OUTPUT output)
if(ipo_supported)
    message("-- IPO is supported!")
    set_property(TARGET tandem PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(WARNING "IPO is not supported: ${output}")
endif()
