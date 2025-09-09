#ifndef FILE_LAPACK_INTERFACE_H
#define FILE_LAPACK_INTERFACE_H

#include <iostream>
#include <string>
#include <vector>
#include <exception>

#include "vector.h"
#include "matrix.h"


#include <complex>

typedef int integer;
typedef integer logical;
typedef float real;
typedef double doublereal;
typedef std::complex<float> singlecomplex;
typedef std::complex<double> doublecomplex;

// Windows SDK defines VOID in the file WinNT.h
#ifndef VOID
typedef void VOID;
#endif
typedef int ftnlen;
typedef int L_fp;  // ?


extern "C" {
#include <clapack.h>
}




namespace ASC_bla
{
  
  // BLAS-1 functions:

  /*
    int daxpy_(integer *n, doublereal *da, doublereal *dx, 
    integer *incx, doublereal *dy, integer *incy);
  */
  // y += alpha x
  template <typename SX, typename SY>
  void AddVectorLapack (double alpha, VectorView<double,SX> x, VectorView<double,SY> y)
  {
    integer n = x.Size();
    integer incx = x.Dist();
    integer incy = y.Dist();
    int err = 
      daxpy_ (&n, &alpha, &x(0),  &incx, &y(0), &incy);
    /* if (err != 0)
      throw std::runtime_error(std::string("daxpy returned errcode "+std::to_string(err)));   */    
  }
  
  
  // BLAS-2 functions:

  // BLAS-3 functions:
  
  // int dgemm_ (char *transa, char *transb, integer *m, integer * n,
  // integer *k, doublereal *alpha, doublereal *a, integer *lda, 
  // doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
  // integer *ldc);

  
  // c = a*b
  template <Ordering OA, Ordering OB>
  void MultMatMatLapack (MatrixView<double, OA> a,
                         MatrixView<double, OB> b,
                         MatrixView<double, Ordering::ColMajor> c)
  {
    char transa_ = (OA == Ordering::ColMajor) ? 'N' : 'T';
    char transb_ = (OB == Ordering::ColMajor) ? 'N' : 'T'; 
  
    integer n = c.Height();
    integer m = c.Width();
    integer k = a.Width();
  
    double alpha = 1.0;
    double beta = 0;
    integer lda = std::max(a.Dist(), 1ul);
    integer ldb = std::max(b.Dist(), 1ul);
    integer ldc = std::max(c.Dist(), 1ul);

    int err =
      dgemm_ (&transa_, &transb_, &n, &m, &k, &alpha, 
              a.Data(), &lda, b.Data(), &ldb, &beta, c.Data(), &ldc);
  }
                       
  template <Ordering OA, Ordering OB>
  void MultMatMatLapack (MatrixView<double, OA> a,
                        MatrixView<double, OB> b,
                        MatrixView<double, Ordering::RowMajor> c)
  {
    MultMatMatLapack (Transpose(b), Transpose(a), Transpose(c));
  }


  

  
  template <Ordering ORD>
  class LapackLU {
    Matrix <double, ORD> a;
    std::vector<integer> ipiv;
  public:
    LapackLU (Matrix<double,ORD> _a)
      : a(std::move(_a)), ipiv(a.Height()) {
      integer m = a.Height();
      if (m == 0) return;
      integer n = a.Width();
      integer lda = a.Dist();
      integer info;
    
      // int dgetrf_(integer *m, integer *n, doublereal *a, 
      //             integer * lda, integer *ipiv, integer *info);
      // std::cout << a;
      dgetrf_(&n, &m, a.Data(), &lda, ipiv.data(), &info);
    }
    
    // b overwritten with A^{-1} b
    void Solve (VectorView<double> b) const {
      char transa =  (ORD == Ordering::ColMajor) ? 'N' : 'T';
      integer n = a.Height();
      integer nrhs = 1;
      integer lda = a.Dist();
      integer ldb = b.Size();
      integer info;

      // int dgetrs_(char *trans, integer *n, integer *nrhs, 
      //             doublereal *a, integer *lda, integer *ipiv,
      //             doublereal *b, integer *ldb, integer *info);

      dgetrs_(&transa, &n, &nrhs, a.Data(), &lda, (integer*)ipiv.data(), b.Data(), &ldb, &info);
    }
  
    Matrix<double,ORD> Inverse() { // rvalue reference return removed
      double hwork;
      integer lwork = -1;
      integer n = a.Height();      
      integer lda = a.Dist();
      integer info;

      // int dgetri_(integer *n, doublereal *a, integer *lda, 
      //             integer *ipiv, doublereal *work, integer *lwork, 
      //             integer *info);

      // query work-size
      dgetri_(&n, &a(0,0), &lda, ipiv.data(), &hwork, &lwork, &info);
      lwork = integer(hwork);
      std::vector<double> work(lwork);
      dgetri_(&n, &a(0,0), &lda, ipiv.data(), &work[0], &lwork, &info);
      return std::move(a);      
    }

    Matrix<double,ORD> LFactor() const { 
      Matrix<double, ORD> l(a.Height(),a.Width());
      l=a;
      for(int i=0;i<std::min(a.Height(),a.Width());i++){
        l(i,i)=1;
      }
      for(int i=0;i<a.Height();i++){
        for(int j=i+1;j<a.Width();j++){
          l(i,j)=0;
        }
      }
      return l;
     }
    Matrix<double,ORD> UFactor() const {  
      Matrix<double, ORD> u(a.Height(),a.Width());
      u=a;
      for(int i=0;i<a.Height();i++){
        for(int j=0;j<i;j++){
          u(i,j)=0;
        }
      }
      return u; 
      }

    std::vector<int> permutations() const{
      int n = std::min(a.Height(),a.Width());
      std::vector<int> per(n);
      for(int i=0;i<n;i++){
        per[i]=i;
      }
      for(int i=0;i<n;i++){
        std::swap(per[i],per[ipiv[i]-1]);
      }
      return per;
    }
    Matrix<double,ORD> PFactor() const { 
      Matrix<double, ORD> p(a.Height(),a.Width());
      std::vector<int> per = permutations();
      for(int i=0;i<std::min(a.Height(),a.Width());i++){
        p(i,per[i])=1;
      }
      return p;
     }
  };
  


  template <Ordering ORD>
  class LapackEigenvalues {
    Matrix <double, ORD> a;
  public:
    LapackEigenvalues (Matrix<double,ORD> _a) : a(_a) {
      if (a.Width() != a.Height()) {
        throw std::invalid_argument("Cannot find eigenvalues of non-quadratic matrix!");
      }
    }

    std::vector<double> SymEigenvalues() {
      char jobz = 'N';
      char uplo = 'U'; // upper triangle
      integer n = a.Width();
      integer lda = a.Dist();
      integer info;
      std::vector<double> eigenvalues(n);

      /* Subroutine
      int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a, 
	      integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	      integer *info);
      */

      // query worksize
      double out_worksize;
      int lwork_query = -1;
      dsyev_(&jobz, &uplo, &n, &a(0,0), &lda, &eigenvalues[0], &out_worksize, &lwork_query, &info);
      int worksize = integer(out_worksize);
      std::vector<double> lwork(worksize);
      
      dsyev_(&jobz, &uplo, &n, &a(0,0), &lda, &eigenvalues[0], &lwork[0], &worksize, &info);

      if (info != 0) {
        throw std::invalid_argument("error calculating eigenvalues!");
      } else {
        return eigenvalues;
      }
    }
  };
  
}


#endif
