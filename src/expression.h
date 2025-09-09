#ifndef FILE_EXPRESSION_H
#define FILE_EXPRESSION_H
#include <exception>
#include <iostream>
#include <utility>
#include <autodiff.h>

namespace ASC_bla {

    template<typename T>
    class VecExpr {
    public:
        auto Upcast() const { return static_cast<const T &> (*this); }

        size_t Size() const { return Upcast().Size(); }

        auto operator()(size_t i) const { return Upcast()(i); }
    };

    template<typename T>
    class MatrixExpr {
    public:
        auto Upcast() const { return static_cast<const T &> (*this); }

        size_t Width() const { return Upcast().Width(); }
        size_t Height() const { return Upcast().Height(); }
        
        auto operator()(size_t row,size_t col) const { return Upcast()(row,col); }
    };

    template<typename T>
    std::ostream& operator<<(std::ostream& oss, const MatrixExpr<T>& m) {            
        for (int i = 0; i < m.Height(); i++) {
            for (int j = 0; j < m.Width(); j++) {
                oss << m(i,j) << " ";
            }
            oss << std::endl;
        }
        return oss;
    }

    template<typename TA, typename TB>
    class SumMatrixExpr : public MatrixExpr<SumMatrixExpr<TA, TB>>{
        TA a_;
        TB b_;
    public:
        SumMatrixExpr(TA a, TB b):a_(a),b_(b){}
        auto operator()(size_t row, size_t col) const { return a_(row,col) + b_(row,col); }
        size_t Width() const { return a_.Width(); }
        size_t Height() const { return a_.Height(); }
    };

    template<typename TA, typename TB>
    class SubstMatrixExpr : public MatrixExpr<SubstMatrixExpr<TA, TB>> {
        TA a_;
        TB b_;
    public:
        SubstMatrixExpr(TA a, TB b):a_(a),b_(b){}
        auto operator()(size_t row, size_t col) const {
            /*std::cout << "start substract" << "\n";
            std::cout << "row:" << row << " col:" << col << "\n"; 
            //std::cout << "a" << "\n" << a_ << "\n";
            std::cout << "b" << "\n" << b_ << "\n";*/

            return a_(row,col) - b_(row,col); 
        }
        size_t Width() const { return a_.Width(); }
        size_t Height() const { return a_.Height(); }
    };

    template<typename TA, typename TB>
    auto operator-(const MatrixExpr<TA> &m, const MatrixExpr<TB> &n) {
        return SubstMatrixExpr(m.Upcast(), n.Upcast());
    }

    template<typename TA, typename SCAL>
    class ScaleMatrixExpr : public MatrixExpr<ScaleMatrixExpr<TA, SCAL>>{
        TA _a;
        SCAL _c;
    public:
        ScaleMatrixExpr(TA a, SCAL c):_a(a),_c{c}{}
        auto operator()(size_t row, size_t col) const { return _a(row,col)*_c; }
        size_t Width() const { return _a.Width(); }
        size_t Height() const { return _a.Height(); }
    };



    template<typename T>
    auto operator*(double scal, const MatrixExpr<T> &m) {
        return ScaleMatrixExpr(m.Upcast(),scal);
    }


    template<size_t D, typename T>
    auto operator*(AutoDiff<D, double> scal, const MatrixExpr<T> &m) {
        return ScaleMatrixExpr(m.Upcast(),scal);
    }

    template<size_t D, typename T>
    auto operator*(const MatrixExpr<T> &m, AutoDiff<D, double> scal) {
        return ScaleMatrixExpr(m.Upcast(),scal);
    }

    template<typename TA, typename TB>
    auto operator+(const MatrixExpr<TA> &m, const MatrixExpr<TB> &n) {
        return SumMatrixExpr(n.Upcast(), m.Upcast());
    }

    template<typename TA, typename TB>
    class MultiplyMatrixExpr : public MatrixExpr<MultiplyMatrixExpr<TA, TB>>{
        TA a_;
        TB b_;
    public:
        MultiplyMatrixExpr(TA a, TB b):a_(a),b_(b){}
        auto operator()(size_t row, size_t col) const {
            if(a_.Width() != b_.Height()){
                throw std::invalid_argument("Matrix dimension must match for multiplication");
            }
            std::remove_cv_t<std::remove_reference_t<decltype(a_(0,0)*b_(0,0))>> sum = a_(0,0)*b_(0,0) - a_(0,0)*b_(0,0); //TODO: Zero element of T
            for (size_t i = 0; i < a_.Width(); i++) {
                sum = sum + (a_(row, i) * b_(i, col));
            }
            return sum;
        }
        size_t Width() const { return b_.Width(); }
        size_t Height() const { return a_.Height(); }
    };

    template<typename TA, typename TB>
    auto operator*(const MatrixExpr<TA> &a, const MatrixExpr<TB> &b) {
        return MultiplyMatrixExpr(a.Upcast(), b.Upcast());
    }

    template<typename TM, typename TV>
    class MatrixMulVecExpr : public VecExpr<MatrixMulVecExpr<TM, TV>> {
        TM _m;
        TV _v;
    public:
        MatrixMulVecExpr(TM m, TV v) : _m(m), _v(v) {}

        auto operator()(size_t row) const {
            if(_m.Width() != _v.Size()){
                throw std::invalid_argument("Matrix/Vector dimension must match for multiplication");
            }
            std::remove_cv_t<std::remove_reference_t<decltype(_m(0,0)*_v(0))>> sum = _v(0)*_m(0,0)-_v(0)*_m(0,0); //TODO: Zero element of ElemT
            for (size_t i = 0; i < _v.Size(); i++) {
                sum += _m(row, i) * _v(i);
            }
            return sum;
        }

        size_t Size() const { return _m.Height(); }
    };

    template<typename TA, typename TB>
    auto operator*(const MatrixExpr<TA> &m, const VecExpr<TB> &v) {
        return MatrixMulVecExpr(m.Upcast(), v.Upcast());
    }

    template<typename TA>
    class TransposeMatExpr : public MatrixExpr<TransposeMatExpr<TA>>{
        TA _a;
    public:
        TransposeMatExpr(TA a):_a(a){}
        auto operator()(size_t row, size_t col) const { return _a(col,row); }
        size_t Width() const { return _a.Height(); }
        size_t Height() const { return _a.Width(); }
    };

  
    class IdMatExpr : public MatrixExpr<IdMatExpr>{
        int dim;
    public:
        IdMatExpr(int _dim):dim(_dim){}
        auto operator()(size_t row, size_t col) const { return (row==col)?1:0; }
        size_t Width() const { return dim; }
        size_t Height() const { return dim; }
    };
    
    template<typename TA, typename TB>
    class SumVecExpr : public VecExpr<SumVecExpr<TA, TB>> {
        TA a_;
        TB b_;
    public:
        SumVecExpr(TA a, TB b) : a_(a), b_(b) {}

        auto operator()(size_t i) const { return a_(i) + b_(i); }

        size_t Size() const { return a_.Size(); }
    };
    template<typename TSCAL, typename TV>
    class ScaleVecExpr : public VecExpr<ScaleVecExpr<TSCAL, TV>> {
        TSCAL scal_;
        TV vec_;
    public:
        ScaleVecExpr(TSCAL scal, TV vec) : scal_(scal), vec_(vec) {}

        auto operator()(size_t i) const { return scal_ * vec_(i); }

        size_t Size() const { return vec_.Size(); }
    };

    template<typename TA, typename TB>
    auto operator+(const VecExpr<TA> &a, const VecExpr<TB> &b) {
        return SumVecExpr(a.Upcast(), b.Upcast());
    }
    template<typename TA, typename TB>
    auto operator-(const VecExpr<TA> &a, const VecExpr<TB> &b) {
        return SumVecExpr(a.Upcast(), ScaleVecExpr<double, TB>(-1, b.Upcast()));
    }
    
    template<typename T>
    auto operator*(double scal, const VecExpr<T> &v) {
        return ScaleVecExpr(scal, v.Upcast());
    }

    template<typename T>
    auto operator*(const VecExpr<T> &v, double scal) {
        return ScaleVecExpr(scal, v.Upcast());
    }

    template<size_t D, typename T>
    auto operator*(AutoDiff<D, double> scal, const VecExpr<T> &v) {
        return ScaleVecExpr(scal, v.Upcast());
    }
    
    template<size_t D, typename T>
    auto operator*(const VecExpr<T> &v, AutoDiff<D, double> scal) {
        return ScaleVecExpr(scal, v.Upcast());
    }


    template<typename T>
    std::ostream &operator<<(std::ostream &ost, const VecExpr<T> &v) {
        if (v.Size() > 0)
            ost << v(0);
        for (size_t i = 1; i < v.Size(); i++)
            ost << ", " << v(i);
        return ost;
    }
}

#endif
