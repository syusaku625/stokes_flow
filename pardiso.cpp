//##################################################################################
//
// FEM solid analysis
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Enginieering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   pardiso_solver.h
 * @author T. Otani
 */

#include "pardiso_solver.h"

using namespace std;

// #################################################################
/**
 * @brief initialize PARDISO class
 * @param [in] DOF  total DOF
 */
void PARDISO_solver::initialize(const int &DOF)
{
  b=(double *)malloc(DOF*sizeof(double));
  x=(double *)malloc(DOF*sizeof(double));
}

void PARDISO_solver::create_csr_matrix(map<pair<int, int>, double> tmp, int numOfNode)
{
    vector<int> column_coo, row_coo;
    vector<double> value_coo;
    for(auto itr = tmp.begin(); itr!=tmp.end(); itr++){
        row_coo.push_back(itr->first.first);
        column_coo.push_back(itr->first.second);
        value_coo.push_back(itr->second);
    }

    nnz = value_coo.size();

    ptr = (MKL_INT *)malloc((numOfNode+1)*sizeof(MKL_INT));
    index = (MKL_INT *)malloc( column_coo.size()*sizeof(MKL_INT) );
    value = (double *)malloc( value_coo.size()*sizeof(double));

    for (int i = 0; i < numOfNode + 1; i++) {
        ptr[i] = 0;
    }

    for (int i = 0; i < value_coo.size(); i++) {
        value[i] = value_coo[i];
        index[i] = column_coo[i];
        ptr[row_coo[i] + 1]++;
    }

    for (int i = 0; i < numOfNode + 1; i++) {
        ptr[i + 1] += ptr[i];
    }
}
// #################################################################
/**
 * @brief set CSR matrix format
 * @param [in] inb          numbers of adjacent nodes of each node
 * @param [in] numOfNode    number of nodes
 * @param [in] dim          total DOF
 */
//void PARDISO_solver::CSR_initialize(const VECTOR2D<int> &inb,const int &numOfNode,const int &dim)
//{
//  ptr = (MKL_INT *)malloc((numOfNode*dim+1)*sizeof(MKL_INT));
//
//  CSR_ptr_initialize(inb,numOfNode,dim);
//  index = (MKL_INT *)malloc( nnz*sizeof(MKL_INT) );
//  value = (double *)malloc( nnz*sizeof(double));
//
//  CSR_index_initialize(inb,numOfNode,dim);
//}

// #################################################################
/**
 * @brief calc L2 norm
 * @param [in] nump DOF
 * @param [in] x    vector
 */
double PARDISO_solver::vector_norm(const int &nump,const double *x)
{
  double norm=0e0;

  #pragma omp parallel for reduction(+:norm)
  for(int i=0;i<nump;i++){
    norm += x[i] * x[i];
  }
  return sqrt(norm);
}

// #################################################################
/**
 * @brief set ptr (CSR matrix) format
 * @param [in] inb          numbers of adjacent nodes of each node
 * @param [in] numOfNode    number of nodes
 * @param [in] dim          total DOF
 */
//void PARDISO_solver::CSR_ptr_initialize(const VECTOR2D<int> &inb,const int &numOfNode,const int &dim)
//{
//  nnz = 0;
//  for(int i=0;i<dim;i++){
//    for(int ic=0;ic<numOfNode;ic++){
//      ptr[ic+i*numOfNode] = nnz;
//      nnz += inb[ic].size()*dim;
//    }
//  }
//  ptr[dim*numOfNode] = nnz;
//}

// #################################################################
/**
 * @brief set index (CSR matrix) format
 * @param [in] inb          numbers of adjacent nodes of each node
 * @param [in] numOfNode    number of nodes
 * @param [in] dim          total DOF
 */
//void PARDISO_solver::CSR_index_initialize(const VECTOR2D<int> &inb,const int &numOfNode,const int &dim)
//{
//  int tmp = 0;
//  for(int dim2=0;dim2<dim;dim2++){
//    for(int ic=0;ic<numOfNode;ic++){
//      for(int k=0;k<dim;k++){
//        for(int i=0;i<inb[ic].size();i++){
//          index[tmp] = inb[ic][i]+k*numOfNode;
//          tmp++;
//        }
//      }
//    }
//  }
//}

// #################################################################
/**
 * @brief set dirichlet boundary condition in CSR matrix
 * @param [in] numOfNode   number of nodes
 * @param [in] ibd         dirichlet boundary mask function
 */
//void PARDISO_solver::set_CSR_dirichlet_boundary_condition1D(const int &numOfBd,ARRAY2D<int> &ibd)
//{
//
//  #pragma omp parallel for
//  for(int it=0;it<numOfBd;it++){
//    
//    int ic = ibd(it,0);
//
//    for(int i=ptr[ic];i<ptr[ic+1];i++){
//      value[i] = 0e0;
//      if(index[i]==ic) value[i] = 1e0;
//    }
//  }
//}

// #################################################################
/**
 * @brief set dirichlet boundary condition in CSR matrix
 * @param [in] numOfNode   number of nodes
 * @param [in] ibd         dirichlet boundary mask function
 */
//void PARDISO_solver::set_CSR_dirichlet_boundary_condition3D(const int &numOfNode,ARRAY2D<int> &ibd)
//{
//
//  #pragma omp parallel for
//  for(int ic=0;ic<numOfNode;ic++){
//
//    if(ibd(ic,0)!=1){
//      for(int i=ptr[ic];i<ptr[ic+1];i++){
//        value[i] = 0e0;
//        if(index[i]==ic) value[i] = 1e0;
//      }
//    }
//    if(ibd(ic,1)!=1){
//      for(int i=ptr[ic+numOfNode];i<ptr[ic+numOfNode+1];i++){
//        value[i] = 0e0;
//        if(index[i]==ic+numOfNode) value[i] = 1e0;
//      }
//    }
//    if(ibd(ic,2)!=1){
//      for(int i=ptr[ic+numOfNode*2];i<ptr[ic+numOfNode*2+1];i++){
//        value[i] = 0e0;
//        if(index[i]==ic+numOfNode*2) value[i] = 1e0;
//      }
//    }
//  }
//}

// #################################################################
/**
 * @brief calc set CSR matrix (Vector)
 * @param [in] K element stiffness matrix
 * @param [in] ie elements
 * @param [in] numOfNode number of nodes
 * @param [in] numOfElm number of elements
 * @param [in] numOfNodeInElm nubmer of node in each elements
 * @param [in] inb nodes around each node
 */
//void PARDISO_solver::set_CSR_value1D(std::vector<ARRAY2D<double>> &K,const std::vector<ElementType> &element,const int &numOfNode,
//                               const int &numOfElm,const VECTOR2D<int> &inb)
//{
//  int tmp1,tmp2,tmp3;
//
//  #pragma omp parallel for
//  for(int ic=0;ic<nnz;ic++) value[ic] = 0e0;
//
//  for(int ielm=0;ielm<numOfElm;ielm++){
//    for(int p=0;p<element[ielm].node.size();p++){
//
//      tmp1 = element[ielm].node[p];
//      for(int q=0;q<element[ielm].node.size();q++){
//
//        tmp2 = element[ielm].node[q];
//
//        for(int i=ptr[tmp1];i<ptr[tmp1+1];i++){
//          if(tmp2==index[i]){
//            value[i]     += K[ielm](p,q);
//            break;
//          }
//        }
//      }
//    }
//  }
//}

// #################################################################
/**
 * @brief calc set CSR matrix (Vector)
 * @param [in] K element stiffness matrix
 * @param [in] ie elements
 * @param [in] numOfNode number of nodes
 * @param [in] numOfElm number of elements
 * @param [in] numOfNodeInElm nubmer of node in each elements
 * @param [in] inb nodes around each node
 */
//void PARDISO_solver::set_CSR_value2D(std::vector<ARRAY4D<double>> &K,const std::vector<ElementType> &element,const int &numOfNode,
//                               const int &numOfElm,const VECTOR2D<int> &inb)
//{
//  int tmp1,tmp2,tmp3;
//
//  #pragma omp parallel for
//  for(int ic=0;ic<nnz;ic++) value[ic] = 0e0;
//
//  for(int ielm=0;ielm<numOfElm;ielm++){
//    for(int p=0;p<element[ielm].node.size();p++){
//
//      tmp1 = element[ielm].node[p];
//      tmp3 = inb[tmp1].size();
//
//      for(int q=0;q<element[ielm].node.size();q++){
//
//        tmp2 = element[ielm].node[q];
//
//        for(int i=ptr[tmp1];i<ptr[tmp1+1];i++){
//          if(tmp2==index[i]){
//            value[i]     += K[ielm](p,q,0,0);
//            value[i+tmp3]  += K[ielm](p,q,0,1);
//            break;
//          }
//        }
//        for(int i=ptr[tmp1+numOfNode];i<ptr[tmp1+numOfNode+1];i++){
//          if(tmp2==index[i]){
//            value[i]     += K[ielm](p,q,1,0);
//            value[i+tmp3]  += K[ielm](p,q,1,1);
//            break;
//          }
//        }
//      }
//    }
//  }
//}

// #################################################################
/**
 * @brief calc set CSR matrix (Vector)
 * @param [in] K element stiffness matrix
 * @param [in] ie elements
 * @param [in] numOfNode number of nodes
 * @param [in] numOfElm number of elements
 * @param [in] numOfNodeInElm nubmer of node in each elements
 * @param [in] inb nodes around each node
 */
//void PARDISO_solver::set_CSR_value3D(std::vector<ARRAY4D<double>> &K,const std::vector<ElementType> &element,const int &numOfNode,
//                               const int &numOfElm,const VECTOR2D<int> &inb)
//{
//  int tmp1,tmp2,tmp3;
//
//  #pragma omp parallel for
//  for(int ic=0;ic<nnz;ic++) value[ic] = 0e0;
//
//  for(int ielm=0;ielm<numOfElm;ielm++){
//    for(int p=0;p<element[ielm].node.size();p++){
//
//      tmp1 = element[ielm].node[p];
//      tmp3 = inb[tmp1].size();
//
//      for(int q=0;q<element[ielm].node.size();q++){
//
//        tmp2 = element[ielm].node[q];
//
//        for(int i=ptr[tmp1];i<ptr[tmp1+1];i++){
//          if(tmp2==index[i]){
//            value[i]     += K[ielm](p,q,0,0);
//            value[i+tmp3]  += K[ielm](p,q,0,1);
//            value[i+tmp3*2] += K[ielm](p,q,0,2);
//            break;
//          }
//        }
//        for(int i=ptr[tmp1+numOfNode];i<ptr[tmp1+numOfNode+1];i++){
//          if(tmp2==index[i]){
//            value[i]     += K[ielm](p,q,1,0);
//            value[i+tmp3]  += K[ielm](p,q,1,1);
//            value[i+tmp3*2] += K[ielm](p,q,1,2);
//            break;
//          }
//        }
//        for(int i=ptr[tmp1+2*numOfNode];i<ptr[tmp1+2*numOfNode+1];i++){
//          if(tmp2==index[i]){
//            value[i]     += K[ielm](p,q,2,0);
//            value[i+tmp3]  += K[ielm](p,q,2,1);
//            value[i+tmp3*2] += K[ielm](p,q,2,2);
//            break;
//          }
//        }
//      }
//    }
//  }
//}


// #################################################################
/**
 * @brief PARDISO MAIn routine
 * @param [in] n          DOF
 * @param [in] numOfOMP   number of OPENMP
 */
void PARDISO_solver::main(MKL_INT n,const int numOfOMP)
{
  for(int ic=0;ic<n+1;ic++) ptr[ic]+=1;
  for(int ic=0;ic<nnz;ic++) index[ic]+=1;

  int mtype = 11;  /*11 Real unsymmetric matrix */

  /* RHS and solution vectors. */
  double *bs=(double *)malloc(n*sizeof(double));
  double res, res0;
  MKL_INT nrhs = 1;     /* Number of right hand sides. */
 
  /* Internal solver memory pointer pt, */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
  /* or void *pt[64] should be OK on both architectures */
  long int *pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  double ddum;          /* Double dummy */
  MKL_INT idum;         /* Integer dummy. */
  char *uplo;
  /* -------------------------------------------------------------------- */
  /* .. Setup Pardiso control parameters. */
  /* -------------------------------------------------------------------- */
  for (int i=0;i<64;i++) iparm[i] = 0;

  iparm[0] = 1;         /* No solver default */
  iparm[1] = 2;         /* Fill-in reordering from METIS */
  //iparm[2] = numOfOMP;      /* number of openMP core in original Pardiso MKL version requires MKL_NUM_THREADS*/
  iparm[3] = 0;         /* No iterative-direct algorithm */
  iparm[4] = 0;         /* No user fill-in reducing permutation */
  iparm[5] = 0;         /* Write solution into x */
  iparm[6] = 0;         /* Not in use */
  iparm[7] = 2;         /* Max numbers of iterative refinement steps */
  iparm[8] = 0;         /* Not in use */
  iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;        /* Conjugate transposed/transpose solve */
  iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  iparm[13] = 0;        /* Output: Number of perturbed pivots */
  iparm[14] = 0;        /* Not in use */
  iparm[15] = 0;        /* Not in use */
  iparm[16] = 0;        /* Not in use */
  iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1;       /* Output: Mflops for LU factorization */
  iparm[19] = 0;        /* Output: Numbers of CG Iterations */
  iparm[26] = 1;

  maxfct = 1;           /* Maximum number of numerical factorizations. */
  mnum = 1;         /* Which factorization to use. */
  msglvl = 0;           /* Print statistical information in file */
  error = 0;            /* Initialize error flag */
  /* -------------------------------------------------------------------- */
  /* .. Initialize the internal solver memory pointer. This is only */
  /* necessary for the FIRST call of the PARDISO solver. */
  /* -------------------------------------------------------------------- */
  for (int i = 0; i < 64; i++ ) pt[i] = 0;
  /* -------------------------------------------------------------------- */
  /* .. Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization. */
  /* -------------------------------------------------------------------- */
  phase = 11;

  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      &n, value, ptr, index, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

  if ( error != 0 )
  {
    printf ("\nERROR during symbolic factorization: %d", error);
    exit (1);
  }
  // printf ("\nReordering completed ... ");
  // printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
  // printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
  /* -------------------------------------------------------------------- */
  /* .. Numerical factorization. */
  /* -------------------------------------------------------------------- */
  phase = 22;
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      &n, value, ptr, index, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  if ( error != 0 )
  {
    printf ("\nERROR during numerical factorization: %d", error);
    exit (2);
  }
  // printf ("\nFactorization completed ... ");
  /* -------------------------------------------------------------------- */
  /* .. Back substitution and iterative refinement. */
  /* -------------------------------------------------------------------- */
  phase = 33;
  //  Loop over 3 solving steps: Ax=b, AHx=b and ATx=b
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
         &n, value, ptr, index, &idum, &nrhs, iparm, &msglvl, b, x, &error);
  /* -------------------------------------------------------------------- */
  /* .. Termination and release of memory. */
  /* -------------------------------------------------------------------- */
  phase = -1;           /* Release internal memory. */
  PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      &n, &ddum, ptr, index, &idum, &nrhs,
      iparm, &msglvl, &ddum, &ddum, &error);

  free(bs);

  for(int ic=0;ic<n+1;ic++) ptr[ic]-=1;
  for(int ic=0;ic<nnz;ic++) index[ic]-=1;
}

