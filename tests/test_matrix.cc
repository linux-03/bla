#include <iostream>
#include "matrix.h"

namespace bla = ASC_bla;



int main()
{
  size_t n = 10;
  bla::Mat<3,3> mat;
  bla::Matrix<double> m(30,30);
  bla::Matrix<double> q(3,3,{0, 0.01, 0.05,
                               7,    0,    0,
                               0,   10, 0.01});
  srand(7893);

  for (int i = 0; i < 30; i++) {
    for (int j = 0; j < 30; j++) {
      m(i,j) = (int)(10*(rand()/(double)RAND_MAX));
      //q(i,j) = (int)(10*(rand()/(double)RAND_MAX));
    }
  }
  /* q(1,0) = 7;
  q(0,2) = 0.05;
  q(0,1) = 0.01;
  q(2,1) = 10;
  q(2, 2)= 0.01; */


  std::cout << m << std::endl << std::endl;
  std::cout << q;

  //bla::MatrixView<double> view = m.Cols(0,2).Rows(1,1);
  //view = 1000;
  
  //m = bla::Matrix<double,bla::Ordering::RowMajor>(10.0 * m);




  bla::Matrix<double> q_inv = bla::inverse(q);

  std::cout << "Inverse\n";
  std::cout << q;
  std::cout << q;
  std::cout << bla::inverse(q);
  std::cout << q*bla::inverse(q);

  q.Rows(1, 2) = {1, 0, 0,
                  0, 1, 0};

  std::cout << std::endl << q << std::endl;

  q = {1, 0, 0,
       0, 1, 0,
       0, 0, 1};

  std::cout << std::endl << q << std::endl;

  /* std::cout << std::endl << std::endl;
  std::cout << "fastMul" << std::endl;
  std::cout << m*q;
  std::cout << std::endl << std::endl;
  std::cout << fastMul(m,q);

  std::cout << std::endl << "difference" << std::endl;
  std::cout << (m*q) + (-1)*fastMul(m,q) << std::endl; */
}