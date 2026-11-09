// Checks for vectors, matrices, expressions and automatic differentiation.
// Run all: ./test_bla        Run one: ./test_bla <name>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

#include "matrix.h"

using namespace ASC_bla;

namespace {

void check(bool ok, const std::string& what) {
  if (!ok) throw std::runtime_error(what);
}

void checkNear(double got, double want, double tol, const std::string& what) {
  if (!(std::fabs(got - want) <= tol))
    throw std::runtime_error(what + ": got " + std::to_string(got) + ", want " + std::to_string(want));
}

template <typename F>
void checkThrows(F f, const std::string& what) {
  bool threw = false;
  try { f(); } catch (const std::exception&) { threw = true; }
  check(threw, what + " should throw");
}

double entry(size_t i, size_t j, double seed) { return std::sin(seed + 1.3 * i + 0.7 * j); }

template <typename M>
void fill(M& m, double seed) {
  for (size_t i = 0; i < m.Height(); i++)
    for (size_t j = 0; j < m.Width(); j++) m(i, j) = entry(i, j, seed);
}

template <typename A, typename B>
double maxDiff(const A& a, const B& b) {
  double d = 0;
  for (size_t i = 0; i < a.Height(); i++)
    for (size_t j = 0; j < a.Width(); j++) d = std::max(d, std::fabs(a(i, j) - b(i, j)));
  return d;
}

Matrix<double> naiveProduct(const MatrixView<double>& a, const MatrixView<double>& b) {
  Matrix<double> c(a.Height(), b.Width());
  for (size_t i = 0; i < a.Height(); i++)
    for (size_t j = 0; j < b.Width(); j++) {
      double s = 0;
      for (size_t k = 0; k < a.Width(); k++) s += a(i, k) * b(k, j);
      c(i, j) = s;
    }
  return c;
}

// --------------------------------------------------------------------------
// Vectors
// --------------------------------------------------------------------------
void test_vector_arithmetic() {
  Vector<double> x(5), y(5);
  for (size_t i = 0; i < 5; i++) { x(i) = i; y(i) = 10; }
  Vector<double> z = x + 3.0 * y - x;
  for (size_t i = 0; i < 5; i++) checkNear(z(i), 30, 0, "x + 3y - x");
  checkNear(x * y, 100, 0, "inner product");
  checkNear(Norm(y), std::sqrt(500.0), 1e-14, "Norm");
  x.Range(1, 3) = 7;
  check(x(0) == 0 && x(1) == 7 && x(2) == 7 && x(3) == 3, "Range assignment");
  x.Slice(0, 2) = -1;
  check(x(0) == -1 && x(1) == 7 && x(2) == -1 && x(4) == -1, "Slice assignment");
}

void test_vector_size_mismatch_throws() {
  Vector<double> a(5), b(3), c(3);
  checkThrows([&] { a = b + c; }, "assigning a 3-vector expression to a 5-vector");
  checkThrows([&] { a.Range(0, 3) = b + a.Range(0, 2); }, "adding vectors of different sizes into a view");
  checkThrows([&] { double d = a * b; (void)d; }, "inner product of different sizes");
  Matrix<double> m(2, 3), n(3, 2);
  checkThrows([&] { auto e = m + n; (void)e; }, "2x3 + 3x2");
  checkThrows([&] { auto e = m - n; (void)e; }, "2x3 - 3x2");
  checkThrows([&] { b.View() = {1.0, 2.0}; }, "initializer list of the wrong length");
}

void test_inner_product_with_infinity() {
  Vector<double> a(2), b(2);
  a(0) = INFINITY; a(1) = 1; b(0) = 1; b(1) = 1;
  check(std::isinf(a * b) && (a * b) > 0, "Inf + 1 should stay Inf");
}

void test_normal_plane_basis() {
  Vector<double> v = {0.3, -1.2, 0.5};
  auto [b1, b2] = normal_plane_basis(v.View());
  checkNear(b1 * v, 0, 1e-14, "b1 orthogonal to v");
  checkNear(b2 * v, 0, 1e-14, "b2 orthogonal to v");
  checkNear(b1 * b2, 0, 1e-14, "b1 orthogonal to b2");
  checkNear(Norm(b1), 1, 1e-14, "b1 unit");
  checkNear(Norm(b2), 1, 1e-14, "b2 unit");
  Vector<double> zero(3);
  checkThrows([&] { normal_plane_basis(zero.View()); }, "normal_plane_basis of the zero vector");
}

// --------------------------------------------------------------------------
// Matrices
// --------------------------------------------------------------------------
void test_matrix_products_match_naive() {
  for (auto [m, n, l] : {std::array<size_t, 3>{1, 1, 1}, {3, 4, 2}, {5, 8, 12}, {7, 3, 13}}) {
    Matrix<double> a(m, n), b(n, l);
    fill(a, 0.1);
    fill(b, 2.3);
    Matrix<double> c = a * b;
    check(maxDiff(c, naiveProduct(a, b)) < 1e-14, "A*B");
    Vector<double> x(n);
    for (size_t i = 0; i < n; i++) x(i) = std::cos(i + 0.5);
    Vector<double> y = a * x;
    for (size_t i = 0; i < m; i++) {
      double s = 0;
      for (size_t k = 0; k < n; k++) s += a(i, k) * x(k);
      checkNear(y(i), s, 1e-14, "A*x");
    }
  }
}

void test_matrix_product_with_infinity() {
  Matrix<double> a(2, 2), id(2, 2);
  a(0, 0) = INFINITY; a(0, 1) = 1; a(1, 0) = 2; a(1, 1) = 3;
  id(0, 0) = id(1, 1) = 1;
  Matrix<double> c = a * id;
  checkNear(c(1, 1), 3, 0, "row without Inf must be unaffected");
  checkNear(c(1, 0), 2, 0, "row without Inf must be unaffected");
  Vector<double> y(2); y(0) = 0; y(1) = 2;
  Vector<double> z = a * y;
  checkNear(z(1), 6, 0, "matrix-vector row without Inf");
}

void test_product_dimension_mismatch_throws() {
  Matrix<double> a(2, 3), b(2, 3);
  checkThrows([&] { auto e = a * b; (void)e; }, "2x3 times 2x3");
  Vector<double> x(2);
  checkThrows([&] { auto e = a * x; (void)e; }, "2x3 times a 2-vector");
}

void test_transpose_of_subviews() {
  Matrix<double> m(3, 4);
  fill(m, 0.4);
  auto cols = m.Cols(1, 2).transpose();
  check(cols.Height() == 2 && cols.Width() == 3, "shape of a transposed column block");
  for (size_t r = 0; r < 2; r++)
    for (size_t c = 0; c < 3; c++) checkNear(cols(r, c), m(c, 1 + r), 0, "transposed column block");
  auto block = m.Block(1, 1, 2, 2).transpose();
  for (size_t r = 0; r < 2; r++)
    for (size_t c = 0; c < 2; c++) checkNear(block(r, c), m(1 + c, 1 + r), 0, "transposed block");
  Matrix<double, Ordering::ColMajor> cm(3, 4);
  fill(cm, 0.4);
  auto rows = cm.Rows(1, 2).transpose();
  for (size_t r = 0; r < 4; r++)
    for (size_t c = 0; c < 2; c++) checkNear(rows(r, c), cm(1 + c, r), 0, "transposed row block of a column-major matrix");
}

void test_multiply_assign() {
  Matrix<double> a(3, 3), b(3, 3);
  fill(a, 1.0);
  fill(b, 2.0);
  Matrix<double> want = naiveProduct(a, b);
  MatrixView<double> view = a;
  view *= b;
  check(maxDiff(a, want) < 1e-14, "A *= B");
  Matrix<double> wide(3, 4);
  checkThrows([&] { MatrixView<double> v = a; v *= wide; }, "3x3 *= 3x4");
}

void test_matrix_assignment() {
  Matrix<double> a(4, 4);
  fill(a, 0.7);
  Matrix<double> copy = a;
  Matrix<double>& alias = a;
  a = alias;
  check(maxDiff(a, copy) == 0, "self-assignment must keep the data");
  Matrix<double> b(2, 2);
  b = a;
  check(b.Height() == 4 && maxDiff(b, a) == 0, "assignment resizes");
  Matrix<double> moved(1, 1);
  moved = std::move(b);
  check(moved.Height() == 4 && maxDiff(moved, copy) == 0, "move assignment");
  Matrix<double> small(3, 3), wide(3, 4);
  checkThrows([&] { small = wide * Matrix<double>(4, 4); }, "assigning a 3x4 expression to a 3x3 matrix");
  checkThrows([&] { Matrix<double> bad(2, 2, {1.0, 2.0, 3.0}); }, "initializer list of the wrong length");
}

void test_fastmul_all_shapes() {
  for (size_t m : {1, 3, 4, 5, 9})
    for (size_t n : {1, 3, 8})
      for (size_t l : {1, 11, 12, 13, 25}) {
        Matrix<double> a(m, n), b(n, l);
        fill(a, 0.3);
        fill(b, 1.9);
        Matrix<double> c = fastMul(a, b);
        check(maxDiff(c, naiveProduct(a, b)) < 1e-13,
              "fastMul " + std::to_string(m) + "x" + std::to_string(n) + "x" + std::to_string(l));
      }
}

void test_products_through_kernel() {
  // Every combination of layouts, transposes and sub-views, against the naive product.
  for (size_t m : {1, 7, 13, 130})
    for (size_t k : {1, 9, 300})
      for (size_t n : {1, 8, 17, 70}) {
        Matrix<double> a(m, k), bt(n, k);
        Matrix<double, Ordering::ColMajor> b(k, n);
        fill(a, 0.2); fill(b, 1.1); fill(bt, 2.3);
        Matrix<double> want = naiveProduct(a, Matrix<double>(b));
        Matrix<double> c(m, n);
        c = a * b;
        check(maxDiff(c, want) < 1e-12, "row-major times column-major");
        Matrix<double, Ordering::ColMajor> cc(m, n);
        cc = a * bt.transpose();
        check(maxDiff(cc, naiveProduct(a, Matrix<double>(bt.transpose()))) < 1e-12, "times a transposed matrix");
        Matrix<double> big(m + 2, n + 3);
        big.Block(1, 2, m, n) = a * b;
        check(maxDiff(big.Block(1, 2, m, n), want) < 1e-12, "product into a sub-view");
        check(big(0, 0) == 0 && big(m + 1, n + 2) == 0, "product into a sub-view leaves the rest alone");
      }
}

void test_nested_and_aliased_products() {
  const size_t n = 40;
  Matrix<double> a(n, n), b(n, n), c(n, n);
  fill(a, 0.1); fill(b, 0.5); fill(c, 0.9);
  Matrix<double> want = naiveProduct(naiveProduct(a, b), c);
  Matrix<double> got = (a * b) * c;
  check(maxDiff(got, want) < 1e-10, "(A*B)*C");
  Matrix<double> got2 = a * (b * c);
  check(maxDiff(got2, want) < 1e-10, "A*(B*C)");
  Matrix<double> scaled = (2.0 * a) * b;
  check(maxDiff(scaled, 2.0 * naiveProduct(a, b)) < 1e-11, "(2A)*B");
  Matrix<double> x = a;
  x = x * b;
  check(maxDiff(x, naiveProduct(a, b)) < 1e-11, "X = X*B");
  Matrix<double> y = b;
  y = a * y;
  check(maxDiff(y, naiveProduct(a, b)) < 1e-11, "Y = A*Y");
}

void test_inverse() {
  Matrix<double> a(5, 5);
  fill(a, 0.9);
  for (size_t i = 0; i < 5; i++) a(i, i) += 3;
  Matrix<double> inv = inverse(a);
  Matrix<double> id = naiveProduct(a, inv);
  for (size_t i = 0; i < 5; i++)
    for (size_t j = 0; j < 5; j++) checkNear(id(i, j), i == j, 1e-13, "A * inverse(A)");
  Matrix<double> singular(2, 2);
  checkThrows([&] { inverse(singular); }, "inverse of a singular matrix");
}

void test_fixed_size_matrix() {
  Mat<2, 3> a;
  check(a(1, 2) == 0, "Mat starts at zero");
  a(1, 2) = 4;
  Mat<2, 3> b = a;
  b(0, 0) = 1;
  check(b(1, 2) == 4 && a(0, 0) == 0, "Mat copies own their data");
  Mat<2, 3> c;
  c = b;
  check(c(1, 2) == 4 && c(0, 0) == 1, "Mat assignment");
}

void test_rotation_matrices() {
  const double deg = 30, s = std::sin(M_PI / 6), c = std::cos(M_PI / 6);
  Matrix<double> want[3] = {Matrix<double>(3, 3, {1, 0, 0, 0, c, -s, 0, s, c}),
                            Matrix<double>(3, 3, {c, 0, s, 0, 1, 0, -s, 0, c}),
                            Matrix<double>(3, 3, {c, -s, 0, s, c, 0, 0, 0, 1})};
  for (int axis = 0; axis < 3; axis++)
    check(maxDiff(makeRotationMatrix3<double>(axis, deg), want[axis]) < 1e-15,
          "rotation about axis " + std::to_string(axis));
}

// --------------------------------------------------------------------------
// Automatic differentiation
// --------------------------------------------------------------------------
template <typename F>
void checkDerivative(F f, double x, const std::string& what) {
  AutoDiff<1> ad(x, 0);
  AutoDiff<1> r = f(ad);
  const double h = 1e-6;
  const double fd = (f(AutoDiff<1>(x + h)).Value() - f(AutoDiff<1>(x - h)).Value()) / (2 * h);
  checkNear(r.DValue(0), fd, 1e-7 * (1 + std::fabs(fd)), what);
}

void test_autodiff_derivatives() {
  // Lambdas return AD explicitly: an expression returned by value would refer to the
  // lambda's own parameter after it is gone.
  using AD = AutoDiff<1>;
  checkDerivative([](AD x) -> AD { return x * x * x; }, 1.7, "x^3");
  checkDerivative([](AD x) -> AD { return 3.0 / x; }, 1.7, "3/x");
  checkDerivative([](AD x) -> AD { return x / (x * x + 1.0); }, 0.4, "x/(x^2+1)");
  checkDerivative([](AD x) -> AD { return 2.0 - x; }, 0.4, "2-x");
  checkDerivative([](AD x) -> AD { return sqrt(x * x + 1.0); }, 0.8, "sqrt of a temporary");
  checkDerivative([](AD x) -> AD { return exp(x); }, 0.8, "exp");
  checkDerivative([](AD x) -> AD { return log(x); }, 0.8, "log");
  checkDerivative([](AD x) -> AD { return sin(x) * cos(x); }, 0.8, "sin*cos");
  checkDerivative([](AD x) -> AD { return tan(x); }, 0.8, "tan");
  checkDerivative([](AD x) -> AD { return pow(x, 2.5); }, 1.3, "pow(x, 2.5)");
  checkDerivative([](AD x) -> AD { return pow(x, 2.0); }, -2.0, "pow(x, 2) at a negative base");
  AD z(0.0, 0);
  AD p = pow(z, 2.0);
  check(p.Value() == 0 && p.DValue(0) == 0, "pow(0, 2)");
}

void test_autodiff_constructors() {
  AutoDiff<3> a(2.0, 0), b(5.0, size_t(2));
  check(a.DValue(0) == 1 && a.DValue(1) == 0 && b.DValue(2) == 1, "seed index");
  AutoDiff<3> zero{};
  check(zero.Value() == 0 && zero.DValue(0) == 0 && zero.DValue(2) == 0, "value-initialized AutoDiff is zero");
}

void test_autodiff_expressions() {
  // Expression templates must give the same results as the operator-by-operator rules.
  using AD = AutoDiff<3>;
  AD x(2.0, 0), y(3.0, 1), z(0.5, 2);
  AD r = x * y + z * x - y / z;
  checkNear(r.Value(), 1.0, 1e-15, "value");
  checkNear(r.DValue(0), 3.5, 1e-15, "d/dx");
  checkNear(r.DValue(1), 0.0, 1e-15, "d/dy");
  checkNear(r.DValue(2), 14.0, 1e-15, "d/dz");
  AD s = 2.0 * x - 1.0 / z + (x - 1.0) * 3.0 + 4.0 - y;
  checkNear(s.Value(), 4 - 2 + 3 + 4 - 3, 1e-15, "scalar mixes, value");
  checkNear(s.DValue(0), 5.0, 1e-15, "scalar mixes, d/dx");
  checkNear(s.DValue(2), 4.0, 1e-15, "scalar mixes, d/dz");
  checkNear(s.DValue(1), -1.0, 1e-15, "scalar mixes, d/dy");
  AD neg = -(x * y);
  checkNear(neg.DValue(0), -3.0, 1e-15, "unary minus");
  // Aliasing: the result may appear on the right-hand side.
  AD a = x;
  a = a * a + a;
  checkNear(a.Value(), 6.0, 1e-15, "a = a*a + a");
  checkNear(a.DValue(0), 5.0, 1e-15, "d(a*a + a)");
  AD b = y;
  b /= x;
  b += x * z;
  checkNear(b.Value(), 2.5, 1e-15, "compound assignment");
  checkNear(b.DValue(0), -0.75 + 0.5, 1e-15, "compound assignment derivative");
  // Functions of expressions.
  AD f = exp(x * z) + log(x + y) + sqrt(x * x + 1.0) + pow(x - y, 2.0);
  const double want = std::exp(1.0) * 0.5 + 1.0 / 5 + 2 / std::sqrt(5.0) + 2 * (2.0 - 3.0);
  checkNear(f.DValue(0), want, 1e-13, "functions of expressions");
  // Temporaries inside an expression stay valid until it is assigned.
  auto make = [](double v) { return AD(v, 1); };
  AD t = make(2.0) * make(3.0) + make(1.0);
  checkNear(t.Value(), 7.0, 1e-15, "temporaries, value");
  checkNear(t.DValue(1), 6.0, 1e-15, "temporaries, derivative");
}

void test_autodiff_nested() {
  using Inner = AutoDiff<3>;
  using AD2 = AutoDiff<2, Inner>;
  AD2 p(Inner(2.0, 0), 0), q(Inner(3.0, 1), 1);
  AD2 r = p * q + 2.0 * p - q / p;
  checkNear(r.Value().Value(), 8.5, 1e-15, "nested value");
  checkNear(r.DValue(0).Value(), 5.75, 1e-15, "nested first derivative");
  checkNear(r.DValue(0).DValue(0), -0.75, 1e-15, "nested second derivative");
}

void test_autodiff_matrix_products() {
  // Products of AutoDiff matrices must differentiate like the scalar case.
  Matrix<AutoDiff<2>> a(2, 2);
  a(0, 0) = AutoDiff<2>(1.5, 0);
  a(0, 1) = AutoDiff<2>(2.0);
  a(1, 0) = AutoDiff<2>(-1.0, 1);
  a(1, 1) = AutoDiff<2>(0.5);
  Matrix<AutoDiff<2>> sq = a * a;
  // (a*a)(0,0) = a00^2 + a01 a10 -> d/da00 = 2 a00 = 3, d/da10 = a01 = 2
  checkNear(sq(0, 0).DValue(0), 3.0, 1e-15, "d(A^2)_00/d a00");
  checkNear(sq(0, 0).DValue(1), 2.0, 1e-15, "d(A^2)_00/d a10");
}

const std::map<std::string, std::function<void()>>& allTests() {
  static const std::map<std::string, std::function<void()>> tests = {
      {"vector_arithmetic", test_vector_arithmetic},
      {"vector_size_mismatch_throws", test_vector_size_mismatch_throws},
      {"inner_product_with_infinity", test_inner_product_with_infinity},
      {"normal_plane_basis", test_normal_plane_basis},
      {"matrix_products_match_naive", test_matrix_products_match_naive},
      {"matrix_product_with_infinity", test_matrix_product_with_infinity},
      {"product_dimension_mismatch_throws", test_product_dimension_mismatch_throws},
      {"transpose_of_subviews", test_transpose_of_subviews},
      {"multiply_assign", test_multiply_assign},
      {"matrix_assignment", test_matrix_assignment},
      {"fastmul_all_shapes", test_fastmul_all_shapes},
      {"inverse", test_inverse},
      {"products_through_kernel", test_products_through_kernel},
      {"nested_and_aliased_products", test_nested_and_aliased_products},
      {"rotation_matrices", test_rotation_matrices},
      {"fixed_size_matrix", test_fixed_size_matrix},
      {"autodiff_derivatives", test_autodiff_derivatives},
      {"autodiff_constructors", test_autodiff_constructors},
      {"autodiff_matrix_products", test_autodiff_matrix_products},
      {"autodiff_expressions", test_autodiff_expressions},
      {"autodiff_nested", test_autodiff_nested},
  };
  return tests;
}

}  // namespace

int main(int argc, char** argv) {
  int failed = 0;
  for (const auto& [name, test] : allTests()) {
    if (argc > 1 && name != argv[1]) continue;
    try {
      test();
      std::printf("PASS %s\n", name.c_str());
    } catch (const std::exception& e) {
      std::printf("FAIL %s: %s\n", name.c_str(), e.what());
      failed++;
    }
  }
  return failed ? 1 : 0;
}
