#!/bin/sh
mkdir build
cd build
cmake -D EIGEN_DIR=/mnt/d/lib/eigen/include/eigen3\
      -D CMAKE_INSTALL_PREFIX=/mnt/d/work/stokes_flow \
      -D TP_DIR=/mnt/d/lib/TextParser \
      ..
make && make install