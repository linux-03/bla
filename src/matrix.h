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

using namespace ASC_HPC;
using namespace std;

namespace ASC_bla {
    enum class Ordering {
        ColMajor, RowMajor
    };

template <typename T=double, Ordering ORD=Ordering::RowMajor>
class MatrixView;

template<typename T=double, Ordering ORD=Ordering::RowMajor>
class Matrix;

template<size_t HEIGHT, size_t WIDTH, typename T=double, Ordering ORD=Ordering::RowMajor>
class Mat;

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
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)+= m2(row, col);
            }
        }
        return *this;
    }

    template<typename TB>
    MatrixView &operator*=(const MatrixExpr<TB> &m2) {
        if(m2.Height() == Height() || m2.Width() == Width()) {
            throw std::invalid_argument("incompatible shapes");
        }
        Matrix<TB> res (Height(), Width());
        res = m2 * *this;

        return (*this)=res;
    }
    MatrixView &operator*=(T s) {
        Matrix<T> res{Height(), Width()};
        res =  s*(*this) ;

        return (*this)=res;
    }
    MatrixView &operator/=(T s) {
        return ((*this)*=(1./s));
    }

    template<typename TB>
    MatrixView &operator-=(const MatrixExpr<TB> &m2) {
        if(_width != m2.Width() || _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
        }
        for (size_t col = 0; col < _width; col++){
            for(size_t row=0; row < _height; row++){
                (*this)(row, col)-= m2(row, col);
            }
        }
        return *this;
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
               return MatrixView<T, ORD==Ordering::RowMajor ? Ordering::ColMajor : Ordering::RowMajor>(_width, _height, _data, ORD == Ordering::RowMajor ? _width : _height);
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
            std::cout << (*this) << std::endl;
            throw std::invalid_argument("Matrix singular");
        }
        RowSwap(column, pivot_i);
    }


    T& operator()(size_t row, size_t column) {
        return const_cast<T&>(std::as_const(*this)(row,column));
    }

    const T& operator()(size_t row, size_t column) const {
        if ((row < 0) || (row >= Height())){
            throw std::invalid_argument(std::string("invalid Matrixview row index ") + std::to_string(row) + " (height: " + std::to_string(Height()) + ")");
        }
        if ((column < 0) || (column >= Width())){
          throw std::invalid_argument(std::string("invalid Matrixview column index ") + std::to_string(column) + " (width: " + std::to_string(Width()) + ")");
        }
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
        : MatrixView<T, ORD> (height, width, new T[list.size()]) {
        // check if list has the right size
        if (list.size() != _height*_width){
            throw std::invalid_argument("initializer list does not have right length for matrix shape");
            return;
        }else{
            // copy list
            for (size_t i = 0; i < list.size(); i++){
                _data[i] = list.begin()[i];
            }
        }
    }

    template<typename TB>
    Matrix &operator=(const MatrixExpr<TB> &m2) {
        if(_width != m2.Width() && _height != m2.Height()){
            throw std::invalid_argument("Matrix dimension must match for copy");
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
        _width = m2._width;
        _height = m2._height;
        delete[] _data;
        _data = new T[_height*_width]();
        _dist = m2._dist;
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
            : MatrixView<T,ORD>(HEIGHT, WIDTH, _local_data.data(), ORD == Ordering::RowMajor ? WIDTH : HEIGHT ) {
                _local_data.fill(0);
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
        if(_width != m2.Width() && _height != m2.Height()){
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

template <size_t SW>
auto InnerProduct4 (size_t n, double * px0, double * px1, double * px2, double * px3, double * py, size_t dy)
{
    SIMD<double,SW> sum0{0.0};
    SIMD<double,SW> sum1{0.0};
    SIMD<double,SW> sum2{0.0};
    SIMD<double,SW> sum3{0.0};
    for (size_t i = 0; i < n; i++)
    {
        // sum += px[i] * SIMD<double,SW>(py+i*dy);
        SIMD<double,SW> yi(py+i*dy);
        sum0 = FMA(SIMD<double,SW>(px0[i]), yi, sum0);
        sum1 = FMA(SIMD<double,SW>(px1[i]), yi, sum1);
        sum2 = FMA(SIMD<double,SW>(px2[i]), yi, sum2);
        sum3 = FMA(SIMD<double,SW>(px3[i]), yi, sum3);
    }
    return std::tuple(sum0, sum1, sum2, sum3);
}

template <size_t SW>
auto InnerProduct (size_t n, double * px, double * py, size_t dy)
{
  SIMD<double,SW> sum{0.0};
  for (size_t i = 0; i < n; i++)
    {
      // sum += px[i] * SIMD<double,SW>(py+i*dy);
      sum = FMA(SIMD<double,SW>(px[i]), SIMD<double,SW>(py+i*dy), sum);
    }
  return sum;
}


inline Matrix<double> fastMul (MatrixView<double, Ordering::RowMajor> a, MatrixView<double, Ordering::RowMajor> b) {
    if(b.Height()!=a.Width())
        throw std::invalid_argument("Matrix dimension must match for multiplication");

    size_t n = a.Width();
    size_t m = a.Height();
    size_t l = b.Width();

    Matrix<double> res(m,l);
    size_t zeilenrest = m % 4;
    size_t spaltenrest = l % 12;

    for (size_t i = 0; i < m / 4; i++) { // Viererpaare von Zeilen
        size_t zeilendist = a.Dist();
        auto zeilenstart = a.Data() + i * 4 * zeilendist;

        for (size_t j = 0; j < l / 12; j++) {
            // inner prod
            size_t spaltenindex = j * 12;
            auto ptr = b.Data() + spaltenindex;

            auto [sum0, sum1, sum2, sum3] = InnerProduct4<12>(
                n,
                zeilenstart,
                zeilenstart + zeilendist,
                zeilenstart + 2 * zeilendist,
                zeilenstart + 3 * zeilendist,
                ptr,
                b.Dist()
            );

            sum0.Store(res.Data() + (i*4)*res.Dist() + spaltenindex);
            sum1.Store(res.Data() + (i*4+1)*res.Dist() + spaltenindex);
            sum2.Store(res.Data() + (i*4+2)*res.Dist() + spaltenindex);
            sum3.Store(res.Data() + (i*4+3)*res.Dist() + spaltenindex);
        }
        size_t lastindex = (l/12)*12;

        for (size_t j = lastindex; j < lastindex + spaltenrest; j++) {
            // rest spalten mit kleinem simd
            auto [sum0, sum1, sum2, sum3] = InnerProduct4<1>(
                n,
                zeilenstart,
                zeilenstart + zeilendist,
                zeilenstart + 2 * zeilendist,
                zeilenstart + 3 * zeilendist,
                b.Data()+j,
                b.Dist()
            );

            sum0.Store(res.Data() + (i*4)*res.Dist() + j);
            sum1.Store(res.Data() + (i*4+1)*res.Dist() + j);
            sum2.Store(res.Data() + (i*4+2)*res.Dist() + j);
            sum3.Store(res.Data() + (i*4+3)*res.Dist() + j);
        }
    }

    size_t lastrow = (n/4)*4;
    // zeilenrest
    for (size_t i = lastrow; i < lastrow + zeilenrest; i++) { // einzelne rest-zeilen
        size_t zeilendist = a.Dist();
        auto zeilenstart = a.Data() + i * zeilendist;

        for (size_t j = 0; j < l / 12; j++) {
            // inner prod
            size_t spaltenindex = j * 12;
            auto ptr = b.Data() + spaltenindex;

            auto sum = InnerProduct<12>(
                n,
                zeilenstart,
                ptr,
                b.Dist()
            );

            sum.Store(res.Data() + i*res.Dist() + spaltenindex);
        }

        size_t lastindex = (l/12)*12;
        for (size_t j = lastindex; j < lastindex + spaltenrest; j++) {
            // rest spalten mit kleinem simd
            auto sum = InnerProduct<1>(
                n,
                zeilenstart,
                b.Data()+j,
                b.Dist()
            );

            sum.Store(res.Data() + i*res.Dist() + j);
        }
    }


    // Matrix<double> res{a.Height(),b.Width()};
    // for(int i=0;i<res.Height();i++){
    //     for(int j=0;j<res.Width();j++){
    //         auto result = SIMD<double,4>(0.);
    //         int rest = b.Height()%4;
    //         for(int l=0; l<b.Height()/4; l++) {
    //             result = FMA(SIMD<double,4>(a.Data()+a.Dist()*i+4*l), SIMD<double,4>(b(0+4*l,j),b(1+4*l,j),b(2+4*l,j),b(3+4*l,j)), result);
    //         }
    //         double sum = HSum(result);
    //         for(int l=rest; l>0; l--) {
    //             int og_index = b.Height()-rest;
    //             sum += a(i,og_index) * b(og_index,j);
    //         }
    //         res(i,j)=sum;
    //     }
    // }

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

// returns rotation matrix for deg degrees around axis: x: 0 y: 1 z: 2
template<typename T>
Matrix<T> makeRotationMatrix3(int axis, T deg){
    if (axis < 0 || axis > 2) throw std::invalid_argument("axis index out of range");

    // convert degree input to radians
    double rad = degToRad<T>(deg);

    // create rotation matrix around x axis
    Matrix<double> R (3, 3, {1, 0,        0,
                                0, std::cos(rad), -std::sin(rad),
                                0, std::sin(rad), std::cos(rad)});

    // correct the rotation axis
    R.ColSwap(0, axis);
    R.RowSwap(0, axis);

    return R;
}

}
#endif
