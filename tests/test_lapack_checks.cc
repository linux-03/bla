// Checks for the LAPACK wrappers. Needs a LAPACK library at link time.
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "lapack_interface.h"

using namespace ASC_bla;

namespace {

int failures = 0;

void report(bool ok, const std::string& what) {
  std::printf("%s %s\n", ok ? "PASS" : "FAIL", what.c_str());
  if (!ok) failures++;
}

template <Ordering ORD>
Matrix<double, ORD> testMatrix(size_t n) {
  Matrix<double, ORD> a(n, n);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) a(i, j) = std::sin(1.0 + 1.3 * i + 0.7 * j) + (i == j ? 4.0 : 0.0);
  return a;
}

template <Ordering ORD>
void checkLU(const char* name) {
  const size_t n = 6;
  Matrix<double, ORD> a = testMatrix<ORD>(n);
  Vector<double> x(n), b(n);
  for (size_t i = 0; i < n; i++) x(i) = std::cos(0.3 * i);
  for (size_t i = 0; i < n; i++) {
    double s = 0;
    for (size_t j = 0; j < n; j++) s += a(i, j) * x(j);
    b(i) = s;
  }
  LapackLU<ORD> lu(a);
  lu.Solve(b);
  double err = 0;
  for (size_t i = 0; i < n; i++) err = std::max(err, std::fabs(b(i) - x(i)));
  report(err < 1e-13, std::string("LU solve, ") + name + " (error " + std::to_string(err) + ")");

  Matrix<double, ORD> inv = LapackLU<ORD>(a).Inverse();
  double ierr = 0;
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      double s = 0;
      for (size_t k = 0; k < n; k++) s += a(i, k) * inv(k, j);
      ierr = std::max(ierr, std::fabs(s - (i == j)));
    }
  report(ierr < 1e-13, std::string("LU inverse, ") + name);

  Matrix<double, ORD> singular(3, 3);
  bool threw = false;
  try { LapackLU<ORD> bad(singular); } catch (const std::domain_error&) { threw = true; }
  report(threw, std::string("singular matrix is reported, ") + name);
}

template <Ordering OA, Ordering OB, Ordering OC>
void checkMultiply(const char* name) {
  const size_t m = 4, k = 5, l = 3;
  Matrix<double, OA> a(m, k);
  Matrix<double, OB> b(k, l);
  Matrix<double, OC> c(m, l);
  for (size_t i = 0; i < m; i++) for (size_t j = 0; j < k; j++) a(i, j) = std::sin(i + 2.0 * j);
  for (size_t i = 0; i < k; i++) for (size_t j = 0; j < l; j++) b(i, j) = std::cos(3.0 * i + j);
  MultMatMatLapack(a.View(), b.View(), c.View());
  double err = 0;
  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < l; j++) {
      double s = 0;
      for (size_t q = 0; q < k; q++) s += a(i, q) * b(q, j);
      err = std::max(err, std::fabs(c(i, j) - s));
    }
  report(err < 1e-13, std::string("dgemm, ") + name);
}

}  // namespace

int main() {
  checkLU<Ordering::RowMajor>("row-major");
  checkLU<Ordering::ColMajor>("column-major");
  checkMultiply<Ordering::ColMajor, Ordering::ColMajor, Ordering::ColMajor>("col * col -> col");
  checkMultiply<Ordering::RowMajor, Ordering::ColMajor, Ordering::ColMajor>("row * col -> col");
  checkMultiply<Ordering::RowMajor, Ordering::RowMajor, Ordering::RowMajor>("row * row -> row");
  checkMultiply<Ordering::ColMajor, Ordering::RowMajor, Ordering::RowMajor>("col * row -> row");

  Matrix<double, Ordering::RowMajor> s(3, 3, {2, 1, 0, 1, 2, 1, 0, 1, 2});
  auto ev = LapackEigenvalues<Ordering::RowMajor>(s).SymEigenvalues();
  const double want[3] = {2 - std::sqrt(2.0), 2, 2 + std::sqrt(2.0)};
  bool ok = ev.size() == 3;
  for (int i = 0; ok && i < 3; i++) ok = std::fabs(ev[i] - want[i]) < 1e-13;
  report(ok, "symmetric eigenvalues");
  return failures ? 1 : 0;
}
