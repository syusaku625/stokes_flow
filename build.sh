#!/bin/sh
mkdir build
cd build
cmake -D EIGEN_DIR=/mnt/c/share/lib/eigen/include/eigen3\
      -D CMAKE_INSTALL_PREFIX=/home/syusaku625/stokes_flow \
      -D TP_DIR=/mnt/c/share/lib/TextParser \
      -D NLopt_DIR=/mnt/c/share/lib/nlopt \
      ..
make && sudo make install