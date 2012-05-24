cmake_minimum_required (VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE true)

find_library(GSL gsl)
message("GSL library is found at " ${GSL})
set (EXTRA_LIBS ${EXTRA_LIBS} ${GSL})

find_library(GSLCBLAS gslcblas)
message("GSLCBLAS library is found at " ${GSLCBLAS})
set (EXTRA_LIBS ${EXTRA_LIBS} ${GSLCBLAS})

find_library(MATH m)
message("MATH library is found at " ${MATH})
set (EXTRA_LIBS ${EXTRA_LIBS} ${MATH})

add_executable (Bessel bessel.cpp)
target_link_libraries (Bessel ${EXTRA_LIBS})

find_library(YAMLCPP yaml-cpp)
if(${YAMLCPP} MATCHES NOTFOUND)
  message("We could not find yaml-cpp: " ${YAMLCPP})
else(${YAMLCPP} MATCHES NOTFOUND)
  message("We use find_library to find yaml-cpp at " ${YAMLCPP})
endif(${YAMLCPP} MATCHES NOTFOUND)


