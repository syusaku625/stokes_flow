#ifndef _PARDISO_SOLVER_H_
#define _PARDISO_SOLVER_H_

//##################################################################################
//
// LIS solver
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   PARDISO_solver.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include <string>
#include <omp.h>
#include <cstdio>
#include<vector>
#include <cstdlib>
#include <cmath>
#include<map>
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include<iostream>

class PARDISO_solver{
 public:
    PARDISO_solver(){
    };
    ~PARDISO_solver(){
      free(ptr);
      free(index);
      free(value);
      free(b);
      free(x);
    };

  //CSR form
  int nnz;
  MKL_INT *ptr,*index;
  double *value,*b,*x;

  void main(MKL_INT n,const int numOfOMP);
  void initialize(const int &DOF);
  //void CSR_initialize(const VECTOR2D<int> &inb,const int &numOfNode,const int &dim);
  void create_csr_matrix(std::map<std::pair<int, int>, double>, int numOfNode);
  double vector_norm(const int &nump,const double *x);

  //void set_CSR_value1D(std::vector<ARRAY2D<double>> &K,const std::vector<ElementType> &element,const int &numOfNode,const int &numOfElm,const VECTOR2D<int> &inb);
  //void set_CSR_value2D(std::vector<ARRAY4D<double>> &K,const std::vector<ElementType> &element,const int &numOfNode,const int &numOfElm,const VECTOR2D<int> &inb);
  //void set_CSR_value3D(std::vector<ARRAY4D<double>> &K,const std::vector<ElementType> &element,const int &numOfNode,const int &numOfElm,const VECTOR2D<int> &inb);
  //void set_CSR_dirichlet_boundary_condition3D(const int &numOfNode,ARRAY2D<int> &ibd);
  // /void set_CSR_dirichlet_boundary_condition1D(const int &numOfBd,ARRAY2D<int> &ibd);

  void main(MKL_INT n,const int numOfOMP,double *b1,double *b2,double *b3);

private:
  //void CSR_ptr_initialize(const VECTOR2D<int> &inb,const int &numOfNode,const int &dim);
  //void CSR_index_initialize(const VECTOR2D<int> &inb,const int &numOfNode,const int &dim);

};

#endif //_PARDISO_SOLVER_H_
