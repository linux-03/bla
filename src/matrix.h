#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>
#include <stdexcept>
#include <vector>
#include <exception>
#include <utility>
#include <initializer_list>
#include <array>
#include <math.h>

#include "expression.h"
#include "vector.h"
#include "simd.h"
#include "gemm.h"
#include <type_traits>

namespace ASC_bla {
    using ASC_HPC::SIMD;
    using ASC_HPC::FMA;
    enum class Ordering {
        ColMajor, RowMajor
    };

template <typename T=double, Ordering ORD=Ordering::RowMajor>
class MatrixView;

template<typename T=double, Ordering ORD=Ordering::RowMajor>
class Matrix;

template<size_t HEIGHT, size_t WIDTH, typename T=double, Ordering ORD=Ordering::RowMajor>
class Mat;

template <typename T>
struct is_double_product : std::false_type {};

template <Ordering OA, Ordering OB>
struct is_double_product<MultiplyMatrixExpr<MatrixView<double,OA>, MatrixView<double,OB>>> : std::true_type {};

template <Ordering OA, Ordering OB>
struct is_double_product<MultiplyMatrixExpr<MatrixView<double,OA>, Matrix<double,OB>>> : std::true_type {};

template <Ordering OA, Ordering OB>
struct is_double_product<MultiplyMatrixExpr<Matrix<double,OA>, MatrixView<double,OB>>> : std::true_type {};

template <Ordering OA, Ordering OB>
struct is_double_product<MultiplyMatrixExpr<Matrix<double,OA>, Matrix<double,OB>>> : std::true_type {};

// c = a * b through the blocked kernel; defined after Matrix.
template <Ordering OC, typename TPROD>
void AssignProduct(MatrixView<double, OC> c, const TPROD & product);

template <typename T, Ordering ORD>
class MatrixView : public MatrixExpr<MatrixView<T,ORD>>{
protected:
    T* _data;
    size_t _width;
    size_t _height;
    size_t _dist;
public:
    MatrixView(size_t height, size_t width, T* data, size_t dist):
        _data{data},_width{width},_height{height},_dist{dist} {};

    MatrixView(size_t height, size_t width, T* data):
        _data{data},_width{width},_height{height},_dist{ORD == Ordering::RowMajor ? width : height} {};

    MatrixView(const MatrixView& m): // copy ctor
        _data{m._data},_width{m._width},_height{m._height},_dist{m._dist} {};

    MatrixView(const MatrixView&& m): // mov ctor
        _data{m._data},_width{m._width},_height{m._height},_dist{m._dist} {};

    MatrixView &operator=(const MatrixView &m2) { //copy assign
        if(_width != m2.Width() || _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)= m2(row, col);
            }
        }
        return *this;
    }

    MatrixView &operator=(MatrixView &&m2) { //mov assign
        if(_width != m2.Width() || _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)= m2(row, col);
            }
        }
        return *this;
    }

    template<typename TB>
    MatrixView &operator=(const MatrixExpr<TB> &m2) {
        if(_width != m2.Width() || _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        if constexpr (std::is_same_v<T, double> && is_double_product<TB>::value) {
            AssignProduct(*this, m2.Upcast());
            return *this;
        }
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)= m2(row, col);
            }
        }
        return *this;
    }

    MatrixView & operator= (std::initializer_list<T> list) {
        if (list.size() != _width*_height){
            throw std::invalid_argument("initializer list does not have right length for matrix shape");
            return *this;
        }
        else{
        for (size_t i = 0; i < _height; i++) {
            for (size_t j = 0; j < _width; j++) {
            if constexpr (ORD == Ordering::RowMajor) {
                _data[_dist * i + j] = list.begin()[_width*i + j];
            } else {
                _data[_dist * j + i] = list.begin()[_height*j + i];
            }
            }
        }
        }
        return *this;
    }

    template<typename TB>
    MatrixView &operator+=(const MatrixExpr<TB> &m2) {
        if(_width != m2.Width() || _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        return *this = (*this) + m2;
    }

    template<typename TB>
    MatrixView &operator*=(const MatrixExpr<TB> &m2) {
        if(m2.Height() != Width() || m2.Width() != Width()) {
            throw std::invalid_argument("incompatible shapes");
        }
        // Evaluate first: writing into *this while reading it would corrupt the result.
        Matrix<T, ORD> product = (*this) * m2;
        return (*this) = product;
    }
    MatrixView &operator*=(T s) {
        return (*this) = s * (*this);
    }
    MatrixView &operator/=(T s) {
        return *this = (1./s) * (*this);
    }

    template<typename TB>
    MatrixView &operator-=(const MatrixExpr<TB> &m2) {
        if(_width != m2.Width() || _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        return *this = (*this) - m2;
    }


    MatrixView &operator=(T scal) {
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)= scal;
            }
        }
        return *this;
    }

    auto Diag() {
        int size = std::min(Width(),Height());
        return VectorView<T, size_t>( size, _dist+1, _data);
    }

    auto View() const { return MatrixView(_height,_width,_data, _dist); }

    size_t Height() const { return _height; }
    size_t Width() const { return _width; }
    size_t Dist() const { return _dist; }
    T* Data(){return _data;}
    const T* Data() const {return _data;}
    MatrixView Cols(size_t first, size_t next) {
        // constructor: size_t height, size_t width, T* data, size_t dist
        return MatrixView(
            this->Height(),
            next,
            &(*this)(0, first),
            this->_dist
        );
    }

    MatrixView Rows(size_t first, size_t next) {
        return MatrixView(
            next,
            this->Width(),
            &(*this)(first, 0),
            this->_dist
        );
    }

    MatrixView Block(size_t first_row, size_t first_col, size_t num_rows, size_t num_cols) {
        return MatrixView(num_rows, num_cols, &(*this)(first_row, first_col), this->_dist);
    }

    VectorView<T, size_t> Row(size_t i) {
        if(ORD == Ordering::RowMajor) {
            return VectorView<T, size_t>(_width, 1, _data+i*_dist);
        }

        return VectorView<T, size_t>(_width, _dist, _data+i);
    }

    void setConstant(T scal) {
        for(size_t i=0; i<Height(); i++) {
            for(size_t j=0; j < Width(); j ++) {
                (*this)(i, j) = scal;
            }
        }
    }

    void diagonal(T scal) {
        for(size_t i = 0; i <Height(); i++) {
            (*this)(i, i) = scal;
        }
    }

    VectorView<T, size_t> Col(size_t j) {
        if(ORD == Ordering::ColMajor) {
            return VectorView<T, size_t>(_height, 1, _data+j*_dist);
        }

        return VectorView<T, size_t>(_height, _dist, _data+j);
    }

    MatrixView<T, ORD==Ordering::RowMajor ? Ordering::ColMajor : Ordering::RowMajor> transpose() {
        // Keep the stride: sub-views (Rows, Cols, Block) have _dist larger than their width.
        return MatrixView<T, ORD==Ordering::RowMajor ? Ordering::ColMajor : Ordering::RowMajor>(_width, _height, _data, _dist);
    };
    void RowSwap(size_t row,size_t dest)  {
        for(size_t i=0; i<Width();i++)
            std::swap((*this)(dest, i),(*this)(row, i));
    }

    void ColSwap(size_t col, size_t dest) {
        for(size_t i=0; i<Height(); i++) {
            std::swap((*this)(i, dest), (*this)(i, col));
        }
    }
    void EnsureNonzero(size_t row, size_t column){
        size_t pivot_i = row;
        T pivot = std::abs((*this)(row,column));
        for(size_t i = row+1; i<Height();i++) {
            if(std::abs((*this)(i,column)) > pivot) {
                pivot = std::abs((*this)(i, column));
                pivot_i = i;
            }
        }
        if(pivot == 0) {
            //std::cout << (*this) << std::endl;
            throw std::invalid_argument("Matrix singular");
        }
        RowSwap(column, pivot_i);
    }


    T& operator()(size_t row, size_t column) {
        return const_cast<T&>(std::as_const(*this)(row,column));
    }

    const T& operator()(size_t row, size_t column) const {
#ifndef NDEBUG
        // Checked only in debug builds: these checks made release code several times slower.
        if (row >= Height()){
            throw std::invalid_argument(std::string("invalid Matrixview row index ") + std::to_string(row) + " (height: " + std::to_string(Height()) + ")");
        }
        if (column >= Width()){
          throw std::invalid_argument(std::string("invalid Matrixview column index ") + std::to_string(column) + " (width: " + std::to_string(Width()) + ")");
        }
#endif
        if constexpr (ORD == Ordering::RowMajor) {
            return _data[row * _dist + column];
        } else {
            return _data[column * _dist + row];
        }
    }
};

template<typename T, Ordering ORD>
class Matrix : public MatrixView<T,ORD> {
    typedef MatrixView<T, ORD> BASE;
    using BASE::_width;
    using BASE::_height;
    using BASE::_data;
    using BASE::_dist;

public:
    Matrix(size_t height, size_t width)
            : MatrixView<T,ORD>(height, width, new T[height*width](), ORD == Ordering::RowMajor ? width : height ) {

            }

    Matrix(const Matrix &m)
            : Matrix(m.Height(),m.Width()) {
        *this = m;
    }

    // initializer list constructor
    Matrix (size_t height, size_t width, std::initializer_list<T> list)
        : MatrixView<T, ORD> (height, width, new T[height*width]()) {
        // check if list has the right size
        if (list.size() != _height*_width){
            delete[] _data;
            throw std::invalid_argument("initializer list does not have right length for matrix shape");
        }else{
            // copy list
            for (size_t i = 0; i < list.size(); i++){
                _data[i] = list.begin()[i];
            }
        }
    }

    template<typename TB>
    Matrix &operator=(const MatrixExpr<TB> &m2) {
        if(_width != m2.Width() || _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        if constexpr (std::is_same_v<T, double> && is_double_product<TB>::value) {
            AssignProduct(MatrixView<T,ORD>(*this), m2.Upcast());
            return *this;
        }
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)= m2(row, col);
            }
        }
        return *this;
    }

    Matrix(Matrix &&m)
            : MatrixView<T,ORD>{0,0, nullptr,0} {
        std::swap(_width, m._width);
        std::swap(_height, m._height);
        std::swap(_data, m._data);
        std::swap(_dist,m._dist);
    }

    template<typename TA>
    Matrix(const MatrixExpr<TA> &m)
            : Matrix(m.Height(),m.Width()) {
        *this = m;
    }


    ~Matrix() {
        delete[] _data;
    }

    Matrix &operator=(const Matrix &m2) {
        if (this == &m2) return *this;
        if (_width*_height != m2._width*m2._height) {
            delete[] _data;
            _data = new T[m2._height*m2._width]();
        }
        _width = m2._width;
        _height = m2._height;
        _dist = m2._dist;
        for (size_t i = 0; i < _width*_height; i++){
            _data[i]=m2._data[i];
        }
        return *this;
    }

    Matrix &operator=(Matrix &&m) noexcept {
        std::swap(_height, m._height);
        std::swap(_width, m._width);
        std::swap(_data, m._data);
        std::swap(_dist, m._dist);
        return *this;
    }

    Matrix & operator= (std::initializer_list<T> list) {
        if (list.size() != _width*_height){
            throw std::invalid_argument("initializer list does not have right length for matrix shape");
            return *this;
        }
        else{
        for (size_t i = 0; i < _height; i++) {
            for (size_t j = 0; j < _width; j++) {
            if constexpr (ORD == Ordering::RowMajor) {
                _data[_dist * i + j] = list.begin()[_width*i + j];
            } else {
                _data[_dist * j + i] = list.begin()[_height*j + i];
            }
            }
        }
        }
        return *this;
    }
};

template<size_t HEIGHT, size_t WIDTH,typename T, Ordering ORD>
class Mat : public MatrixView<T,ORD> {
    typedef MatrixView<T, ORD> BASE;
    using BASE::_width;
    using BASE::_height;
    using BASE::_data;
    using BASE::_dist;
    std::array<T, WIDTH*HEIGHT> _local_data;

public:
    Mat()
            : MatrixView<T,ORD>(HEIGHT, WIDTH, nullptr, ORD == Ordering::RowMajor ? WIDTH : HEIGHT ) {
                // Point the view at the array only once the array exists.
                _data = _local_data.data();
                _local_data.fill(T{});
            }

    Mat(const Mat &m)
            : Mat() {
        *this = m;
    }

    // initializer list constructor
    Mat (std::initializer_list<T> list)
        : Mat() {
        // check if list has the right size
        if (list.size() != _height*_width){
            throw std::invalid_argument("initializer list does not have right length for matrix shape");
            return;
        }else{
            // copy list
            for (size_t i = 0; i < list.size(); i++){
                _local_data[i] = list.begin()[i];
            }
        }
    }

    template<typename TB>
    Mat &operator=(const MatrixExpr<TB> &m2) {
        if(_width != m2.Width() || _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)= m2(row, col);
            }
        }
        return *this;
    }

    template<typename TA>
    Mat(const MatrixExpr<TA> &m)
            : Mat() {
        *this = m;
    }


    Mat &operator=(const Mat &m2) {
        for (size_t i = 0; i < _width*_height; i++){
            _data[i]=m2._data[i];
        }
        return *this;
    }

    // Matrix &operator=(Matrix &&m) {
    //     std::swap(_height, m._height);
    //     std::swap(_width, m._width);
    //     std::swap(_data, m._data);
    //     std::swap(_dist, m._dist);
    //     return *this;
    // }

    Mat & operator= (std::initializer_list<T> list) {
        if (list.size() != _width*_height){
            throw std::invalid_argument("initializer list does not have right length for matrix shape");
            return *this;
        }
        else{
        for (size_t i = 0; i < _height; i++) {
            for (size_t j = 0; j < _width; j++) {
            if constexpr (ORD == Ordering::RowMajor) {
                _data[_dist * i + j] = list.begin()[_width*i + j];
            } else {
                _data[_dist * j + i] = list.begin()[_height*j + i];
            }
            }
        }
        }
        return *this;
    }
};

template<typename T, Ordering ORD>
MatrixView<T, ORD==Ordering::RowMajor ? Ordering::ColMajor : Ordering::RowMajor> transpose(MatrixView<T, ORD> &m) {
    return m.transpose();
}

template <typename T>
auto Evaluated(const MatrixExpr<T> & m) {
    using Elem = std::remove_cv_t<std::remove_reference_t<decltype(m.Upcast()(0,0))>>;
    if constexpr (std::is_base_of_v<MatrixView<Elem, Ordering::RowMajor>, T> ||
                  std::is_base_of_v<MatrixView<Elem, Ordering::ColMajor>, T>)
        return m.Upcast().View();
    else
        return Matrix<Elem>(m.Upcast());
}

namespace detail {
    // Strides of element (i, j): data[i * rs + j * cs].
    template <typename T, Ordering ORD>
    std::pair<ptrdiff_t, ptrdiff_t> Strides(const MatrixView<T, ORD> & m) {
        if constexpr (ORD == Ordering::RowMajor) return { ptrdiff_t(m.Dist()), 1 };
        else return { 1, ptrdiff_t(m.Dist()) };
    }

    // Address range of a view, for the overlap test.
    template <typename T, Ordering ORD>
    std::pair<const T*, const T*> Extent(const MatrixView<T, ORD> & m) {
        if (m.Height() == 0 || m.Width() == 0) return { m.Data(), m.Data() };
        const T * first = m.Data();
        return { first, &m(m.Height()-1, m.Width()-1) + 1 };
    }

    template <typename T, Ordering O1, Ordering O2>
    bool Overlap(const MatrixView<T, O1> & x, const MatrixView<T, O2> & y) {
        auto [x0, x1] = Extent(x);
        auto [y0, y1] = Extent(y);
        return x0 < y1 && y0 < x1;
    }
}

template <Ordering OC, typename TPROD>
void AssignProduct(MatrixView<double, OC> c, const TPROD & product) {
    auto a = product.A().View();
    auto b = product.B().View();
    const size_t m = c.Height(), n = c.Width(), k = a.Width();
    const bool overlap = detail::Overlap(c, a) || detail::Overlap(c, b);
    if (m * n * k <= 64 && !overlap) {
        // Tiny products (3x3 and similar): a direct loop beats the kernel's setup.
        auto [rsa, csa] = detail::Strides(a);
        auto [rsb, csb] = detail::Strides(b);
        auto [rsc, csc] = detail::Strides(c);
        const double * pa = a.Data();
        const double * pb = b.Data();
        double * pc = c.Data();
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++) {
                double sum = 0;
                for (size_t q = 0; q < k; q++) sum += pa[i * rsa + q * csa] * pb[q * rsb + j * csb];
                pc[i * rsc + j * csc] = sum;
            }
        return;
    }
    if (overlap) {
        // c = c * b and similar: compute into fresh storage first.
        Matrix<double, OC> tmp(c.Height(), c.Width());
        AssignProduct(MatrixView<double, OC>(tmp), product);
        c = tmp;
        return;
    }
    auto [rsa, csa] = detail::Strides(a);
    auto [rsb, csb] = detail::Strides(b);
    auto [rsc, csc] = detail::Strides(c);
    Gemm(c.Height(), c.Width(), a.Width(), 1.0, a.Data(), rsa, csa, b.Data(), rsb, csb, 0.0, c.Data(), rsc, csc);
}

// Kept for compatibility: the blocked kernel now handles every shape and layout.
inline Matrix<double> fastMul (MatrixView<double, Ordering::RowMajor> a, MatrixView<double, Ordering::RowMajor> b) {
    if(b.Height()!=a.Width())
        throw std::invalid_argument("Matrix dimension must match for multiplication");
    Matrix<double> res(a.Height(), b.Width());
    res = a * b;
    return res;
}

template <typename T>
Vector<T> AsVector (MatrixView<T> mat)  {
    VectorView<T> res(mat.Height()*mat.Width(), mat.Data());
    return res;
}

template <typename T>
auto AsMatrix (VectorView<T> vec, size_t h, size_t w) {
    return MatrixView<T> (h,w,vec.Data());
}

template <typename T>
Matrix<T> inverse(const MatrixView<T>& m ){
    if(m.Height() != m.Width())
        throw std::invalid_argument("Only square matrixes allowed");
    Matrix<T> work(m.Height(),m.Height()*2);
    for(size_t row=0; row < work.Height();row++){
        for(size_t column=0; column < work.Height();column++){
            work(row,column)=m(row,column);
            work(row,column+m.Height()) = (row==column) ? 1 : 0;
        }
    }
    for(size_t column=0; column < work.Height();column++){
        work.EnsureNonzero(column,column);

        work.Rows(column,1)=1.0/(work(column,column))*work.Rows(column,1);
        //work.RowMultiply(column,1/(work(column,column)));

        for(size_t row=0; row < work.Height();row++){
            if(row==column) continue;
            work.Rows(row, 1) = work.Rows(row, 1) + (-work(row,column))*work.Rows(column,1);
            //work = work.RowMulAdd(column,row,-work(row,column));

        }
    }
    Matrix<T> result = work.Cols(m.Width(),m.Width());
    return result;
}

template<typename T>
Matrix<T> Diagonal(size_t height, T el) {
    Matrix<T> eye(height, height);

    eye.Diag() = el;

    return eye;
}

template<typename T>
T degToRad(T deg){
  return M_PI/180. * deg;
}

template<typename T>
T radToDeg(T deg){
  return 180./M_PI * deg;
}

// Rotation by deg degrees (counter-clockwise, right-hand rule) about axis x: 0, y: 1, z: 2.
template<typename T>
Matrix<T> makeRotationMatrix3(int axis, T deg){
    if (axis < 0 || axis > 2) throw std::invalid_argument("axis index out of range");
    const T rad = degToRad<T>(deg);
    const T c = std::cos(rad), s = std::sin(rad);
    // The other two axes, in cyclic order, span the plane of rotation.
    const int i = (axis + 1) % 3, j = (axis + 2) % 3;
    Matrix<T> R(3, 3);
    R(axis, axis) = 1;
    R(i, i) = c;  R(i, j) = -s;
    R(j, i) = s;  R(j, j) = c;
    return R;
}

}
#endif
