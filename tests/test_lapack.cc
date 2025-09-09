#include <iostream>

#include <vector.h>
#include <lapack_interface.h>
#include <chrono>
#include <cmath>
#include <stdlib.h>

//#define RUN_TIMINGS

using namespace ASC_bla;
using namespace std;


int main()
{
  Vector<double> x(5);
  Vector<double> y(5);

  Matrix<double> q(3,3);
  Matrix<double> p(3,3);
  Matrix<double> c(3,3);
  Matrix<double> a(5,5);

  q(0,0)=1;
  q(0,1)=2;
  q(0,2)=3;
  q(1,1)=4;
  q(1,2)=5;
  q(2,2)=6;
  /* Matrix q
  1 2 3
  0 4 5
  0 0 6
  */
  
  p(0,0)=1;
  p(0,1)=2;
  p(0,2)=3;
  p(1,1)=4;
  p(1,2)=5;
  p(2,2)=6;
  p(2,0)=0;
  /* Matrix p
  1 2 3
  0 4 5
  0 0 6
  */

  for (int i = 0; i < a.Height(); i++) {
    for (int j = 0; j <= i; j++) {
      a(i,j) = 3*i + j/2;
      a(j,i) = 3*i + j/2;
    }
  }
  /* Matrix a 
  0  3  6  9  12
  3  3  6  9  12
  6  6  7  10 13
  9  9  10 10 13
  12 12 13 13 14
  */

  cout << a;

  for (int i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 2;
    }

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  
  AddVectorLapack (2,x,y);
  cout << "y+2*x = " << y << endl;

  MultMatMatLapack(p,q,c);
  std::cout << c;

  #ifdef RUN_TIMINGS
  cout << "lapack" << std::endl;
    srand(0);
  for(int l=1; l<=3;l++){
    int n = pow(10,l);
    Matrix<double> m1(n,n);
    Matrix<double> m2(n,n);
    Matrix<double> m3(n,n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        m1(i,j)=(double)(rand())/RAND_MAX;
        m2(i,j)=(double)(rand())/RAND_MAX;
        m3(i,j)=0;
      }
    }
    size_t flops = n*n*n;
    size_t runs = size_t (1e9 / flops) + 1;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < runs; i++){
      // lapack
      MultMatMatLapack(m1,m2,m3);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end-start).count();
          
    cout << "n = " << n << ", time = " << time << " s, GFlops = " 
        << (flops*runs)/time*1e-9 << endl;
  }

  cout << "fastMul" << std::endl;
  srand(0);
  for(int l=1; l<=3;l++){
    int n = pow(10,l);
    Matrix<double> m1(n,n);
    Matrix<double> m2(n,n);
    Matrix<double> m3(n,n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        m1(i,j)=(double)(rand())/RAND_MAX;
        m2(i,j)=(double)(rand())/RAND_MAX;
        m3(i,j)=0;
      }
    }
    size_t flops = n*n*n;
    size_t runs = size_t (1e9 / flops) + 1;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < runs; i++){
      // fastmul
      m3 = fastMul(m1, m2);
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end-start).count();
          
    cout << "n = " << n << ", time = " << time << " s, GFlops = " 
        << (flops*runs)/time*1e-9 << endl;
  }

  cout << "standard" << std::endl;
  for(int l=1; l<=3;l++){
    int n = pow(10,l);
    Matrix<double> m1(n,n);
    Matrix<double> m2(n,n);
    Matrix<double> m3(n,n);
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
        m1(i,j)=(double)(rand())/RAND_MAX;
        m2(i,j)=(double)(rand())/RAND_MAX;
        m3(i,j)=0;
      }
    }
    size_t flops = n*n*n;
    size_t runs = size_t (1e9 / flops) + 1;

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < runs; i++){
      m3=m1*m2;
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end-start).count();
          
    cerr << m3;

    cout << "n = " << n << ", time = " << time << " s, GFlops = " 
        << (flops*runs)/time*1e-9 << endl;
    }

  #endif

  std::vector<double> eigenvalues = LapackEigenvalues(a).SymEigenvalues();
  cout << endl;
  for (int i = 0; i < eigenvalues.size(); i++) {
    cout << eigenvalues[i];
    cout << ",";


  }
  std::cout << p;
  p=p.Transpose();
  std::cout << p.Transpose();
  std::cout << "L:\n";
  Matrix<double> p2 = q.Transpose();
  std::cout << p2;
  LapackLU<Ordering::RowMajor> lu(p2);
  auto l = lu.LFactor();
  auto u = lu.UFactor();
  auto per = lu.PFactor();
  std::cout<<l;
  std::cout << "U:\n";
  std::cout<< u;
  std::cout << "P:\n";
  std::cout<< per;
  std::cout << "P*L*U:\n";
  std::cout<< l*u;
}
