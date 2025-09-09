#include <iostream>

#include <vector.h>
#include "matrix.h"

namespace bla = ASC_bla;


int main()
{
  // size_t n = 5;
  // bla::Vector<double> x(n), y(n), d(2);
  // bla::Matrix<double> a(2,2), b(2,2);
  // a(0,0)=1;
  // a(0,1)=2;
  // a(1,0)=3;
  // a(1,1)=4;

  // b(0,0)=4;
  // b(0,1)=3;
  // b(1,0)=2;
  // b(1,1)=1;

  // d(0)=3;
  // d(1)=1;

  // for (size_t i = 0; i < x.Size(); i++)
  //   {
  //     x(i) = i;
  //     y(i) = 10;
  //   }

  // a.Dump();
  // printf("\n");
  // b.Dump();
  // printf("\n");
  // (a*b).Dump();
  // bla::Vector<double> z = x+y;

  // d.Dump();
  
  // std::cout << "x+y = " << z << std::endl;
  // (a*d).Dump();
  // std::cout << "\n";
  // (a*3.0).Dump();
  // a.Inverse();
  // b.Inverse();

  size_t n = 10;
  bla::Vector<double> x(n), y(n);

  for (size_t i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 10;
    }

  bla::Vector<double> z = x+y;
  
  std::cout << "x+y = " << z << std::endl;
  
  std::cout << "type of (x+3*y) is  " << typeid(x+3*y).name() << std::endl;

  std::cout << "x+3*y = " << x+3*y << std::endl;

  std::cout << "sizeof(x+3*y) = " << sizeof(x+3*y) << std::endl;
  std::cout <<x<<"\n";
  x.Range(2,9) = 3;
  std::cout <<x<<"\n";
  x.Slice(1,5) = 10;
  std::cout << "x = " << x << std::endl;
  x.Range(1,5) = {11, 12, 13, 14};
  std::cout << "x = " << x << std::endl;

  bla::Vector<double> w({1, 2, 3, 4, 5})
  std::cout << "x = " << x << std::endl;

}


