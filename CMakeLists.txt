###################################################################################
#
# Vessel tree construction
#
# Copyright (c) 2020 Biomechanics Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

cmake_minimum_required(VERSION 3.16)

option(enable_OPENMP "Enable OpenMP" ON)
option(enable_GLOG "Enable GLOG" OFF)

set(CMAKE_C_COMPILER "icc")
set(CMAKE_CXX_COMPILER "icpc")
set(CMAKE_MPI_C_COMPILER "mpiicc")
set(CMAKE_MPI_CXX_COMPILER "mpiicpc")

set(base_compile_options "-Wall -Wextra -O3")
set(base_compile_options ${base_compile_options} " -qopenmp")
string(REPLACE " " " " compile_options ${base_compile_options})
#-----------------------

project(fem_homozinied_diffusion C CXX)
set(CMAKE_CXX_STANDARD 11)

#change default install prefix
set(INSTALL_DIR "/home/syusaku625/bayesopt_test")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH "..." FORCE)
endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

# TextParser
set(TP_DIR ${TP_DIR} CACHE PATH "...")
find_package(TP REQUIRED)
include_directories(${TP_INCLUDE_DIR})

# MKL
set(MKL_ARCH "intel64")
set(MKL_LINK "static")
set(MKL_INTERFACE "lp64")
set(MKL_THREADING "intel_thread")
set(MKL_MPI "intelmpi")
find_package(MKL REQUIRED C CXX)
include_directories(${MKL_INCLUDE_DIR})
MESSAGE(${MKL_LIBRARIES})

add_subdirectory(lib_src)

set(main_target "Main")

add_executable(${main_target} main.cpp)
target_sources(${main_target}
  PRIVATE
    main.cpp
)

target_compile_options(${main_target} PUBLIC ${compile_options_tmp})

target_link_libraries(${main_target} PRIVATE
  STOKES
  ${TP_LIBRARY}
  ${MKL_LIBRARIES}
  iomp5
)

install(TARGETS ${main_target} DESTINATION ${CMAKE_INSTALL_PREFIX})