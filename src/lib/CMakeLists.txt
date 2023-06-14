###################################################################################
#
# Coherent point drift
#
# Copyright (c) 2020 Biomechanics Lab.,
#                    Department of Mechanical Science and Bioengineering,
#                    Graduate School of Engineering Science,
#                    Osaka University.
# All rights reserved.
#
###################################################################################

set(compile_options_tmp ${base_compile_options})
string(REPLACE " " ";" compile_options_tmp ${compile_options_tmp})

add_library(STOKES STATIC)
target_sources(STOKES
  PRIVATE
    pardiso.cpp
    domain.cpp
    gauss.cpp
    fileIO.cpp
    parameter.cpp
    stokes_main.cpp
    topology_stokes.cpp
    MMA.cpp
)

target_include_directories(STOKES PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

target_compile_options(STOKES PUBLIC ${compile_options_tmp})