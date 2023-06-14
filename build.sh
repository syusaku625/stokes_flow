#!/bin/sh
mkdir build
cd build
cmake -D EIGEN_DIR=/mnt/d/lib/eigen/include/eigen3\
      -D CMAKE_INSTALL_PREFIX=/mnt/d/work/stokes_flow/bin \
      -D TP_DIR=/mnt/d/lib/TextParser \
      -D NLopt_DIR=/mnt/d/lib/nlopt \
      ..
make && make install