#ifndef FILE_AUTODIFF_CC
#define FILE_AUTODIFF_CC

// Forward-mode automatic differentiation.
//
// AutoDiff<D, SCAL> holds a value and D partial derivatives. For SCAL = double,
// arithmetic builds expression templates: x*y + z*x - y/z is a small tree of
// nodes that is evaluated in a single loop over the D derivatives when it is
// assigned, instead of creating a full temporary per operation. Nodes cache
// their values when built, and keep AutoDiff operands by reference only when
// those are named objects (temporaries are copied), so expressions stay valid
// until the end of the statement that creates them. Store results in an
// AutoDiff, not in `auto`.
//
// Nested types such as AutoDiff<12, AutoDiff<30>> use direct (eager)
// operators: their values are themselves AutoDiff objects.

#include <cmath>
#include <iostream>
#include <type_traits>
#include <utility>

// Marks a loop whose iterations are independent, so it can be vectorized.
// Left out for WebAssembly, which is built without SIMD and would only warn.
#if defined(__EMSCRIPTEN__)
#define ASC_BLA_AD_INDEPENDENT_ITERATIONS
#elif defined(__clang__)
#define ASC_BLA_AD_INDEPENDENT_ITERATIONS _Pragma("clang loop vectorize(assume_safety)")
#elif defined(__GNUC__)
#define ASC_BLA_AD_INDEPENDENT_ITERATIONS _Pragma("GCC ivdep")
#else
#define ASC_BLA_AD_INDEPENDENT_ITERATIONS
#endif

namespace ASC_bla
{
    template <size_t D, typename SCAL = double> class AutoDiff;

    namespace ad_detail {

        // Marks expression nodes; AutoDiff<D, double> also counts as an expression.
        struct ExprTag {};
        // Base of nested AutoDiff types, which are not expressions.
        struct NoTag {};

        template <typename T>
        constexpr bool is_expr = std::is_base_of_v<ExprTag, std::decay_t<T>>;

        template <typename T>
        struct is_leaf : std::false_type {};
        template <size_t D>
        struct is_leaf<AutoDiff<D, double>> : std::true_type {};

        template <typename T>
        constexpr bool is_node = is_expr<T> && !is_leaf<std::decay_t<T>>::value;

        template <typename T>
        constexpr bool is_scalar = std::is_arithmetic_v<std::decay_t<T>>;

        // Named AutoDiff objects are held by reference; everything else by value.
        template <typename T>
        using Hold = std::conditional_t<std::is_lvalue_reference_v<T> && is_leaf<std::decay_t<T>>::value,
                                        const std::decay_t<T> &, std::decay_t<T>>;

        template <typename A, typename B>
        constexpr size_t SameDim() {
            static_assert(std::decay_t<A>::DIM == std::decay_t<B>::DIM, "AutoDiff dimensions must match");
            return std::decay_t<A>::DIM;
        }

        template <typename A, typename B>
        class Sum : public ExprTag {
            A a; B b; double v;
        public:
            static constexpr size_t DIM = SameDim<A, B>();
            template <typename X, typename Y>
            Sum(X && x, Y && y) : a(std::forward<X>(x)), b(std::forward<Y>(y)), v(a.Value() + b.Value()) {}
            double Value() const { return v; }
            double DValue(size_t i) const { return a.DValue(i) + b.DValue(i); }
        };

        template <typename A, typename B>
        class Diff : public ExprTag {
            A a; B b; double v;
        public:
            static constexpr size_t DIM = SameDim<A, B>();
            template <typename X, typename Y>
            Diff(X && x, Y && y) : a(std::forward<X>(x)), b(std::forward<Y>(y)), v(a.Value() - b.Value()) {}
            double Value() const { return v; }
            double DValue(size_t i) const { return a.DValue(i) - b.DValue(i); }
        };

        template <typename A, typename B>
        class Prod : public ExprTag {
            A a; B b; double av, bv;
        public:
            static constexpr size_t DIM = SameDim<A, B>();
            template <typename X, typename Y>
            Prod(X && x, Y && y) : a(std::forward<X>(x)), b(std::forward<Y>(y)), av(a.Value()), bv(b.Value()) {}
            double Value() const { return av * bv; }
            double DValue(size_t i) const { return a.DValue(i) * bv + av * b.DValue(i); }
        };

        template <typename A, typename B>
        class Quot : public ExprTag {
            A a; B b; double q, inv;
        public:
            static constexpr size_t DIM = SameDim<A, B>();
            template <typename X, typename Y>
            Quot(X && x, Y && y) : a(std::forward<X>(x)), b(std::forward<Y>(y)),
                                   q(a.Value() / b.Value()), inv(1.0 / b.Value()) {}
            double Value() const { return q; }
            // Quotient rule in one pass: (a/b)' = (a' - (a/b) b') / b.
            double DValue(size_t i) const { return (a.DValue(i) - q * b.DValue(i)) * inv; }
        };

        // c * a, and a / c as (1/c) * a.
        template <typename A>
        class Scale : public ExprTag {
            A a; double c, v;
        public:
            static constexpr size_t DIM = std::decay_t<A>::DIM;
            template <typename X>
            Scale(X && x, double cc) : a(std::forward<X>(x)), c(cc), v(cc * a.Value()) {}
            double Value() const { return v; }
            double DValue(size_t i) const { return c * a.DValue(i); }
        };

        // s * a + c, which covers a + c, c + a, a - c, c - a and -a.
        template <typename A>
        class Affine : public ExprTag {
            A a; double s, v;
        public:
            static constexpr size_t DIM = std::decay_t<A>::DIM;
            template <typename X>
            Affine(X && x, double ss, double c) : a(std::forward<X>(x)), s(ss), v(ss * a.Value() + c) {}
            double Value() const { return v; }
            double DValue(size_t i) const { return s * a.DValue(i); }
        };

        // c / a.
        template <typename A>
        class Recip : public ExprTag {
            A a; double v, f;
        public:
            static constexpr size_t DIM = std::decay_t<A>::DIM;
            template <typename X>
            Recip(X && x, double c) : a(std::forward<X>(x)) {
                const double inv = 1.0 / a.Value();
                v = c * inv;
                f = -v * inv;
            }
            double Value() const { return v; }
            double DValue(size_t i) const { return f * a.DValue(i); }
        };

    } // namespace ad_detail


    // AutoDiff<D, double> for an expression node, the type itself otherwise.
    // Generic code uses this to declare variables that hold a result.
    template <typename T, bool = ad_detail::is_node<T>>
    struct ad_value { using type = std::remove_cv_t<std::remove_reference_t<T>>; };
    template <typename T>
    struct ad_value<T, true> { using type = AutoDiff<std::decay_t<T>::DIM, double>; };
    template <typename T>
    using ad_value_t = typename ad_value<T>::type;


    template <size_t D, typename SCAL>
    class AutoDiff : public std::conditional_t<std::is_same_v<SCAL, double>, ad_detail::ExprTag, ad_detail::NoTag>
    {
        SCAL val;
        SCAL dval[D?D:1];

        template <typename E>
        static constexpr bool is_node = ad_detail::is_node<E>;

    public:
        static constexpr size_t DIM = D;
        typedef AutoDiff<D, SCAL> TELEM;
        typedef SCAL TSCAL;

        /// Uninitialized; AutoDiff{} (value initialization) is zero.
        AutoDiff () = default;
        AutoDiff (const AutoDiff &) = default;
        AutoDiff & operator= (const AutoDiff &) = default;

        /// constant value
        AutoDiff (SCAL aval) : val(aval)
        {
            for (size_t i = 0; i < D; i++) dval[i] = SCAL(0);
        }

        /// constant value of another arithmetic or convertible type
        template <typename T, typename std::enable_if<!ad_detail::is_expr<T> && std::is_convertible_v<T, SCAL>, int>::type = 0>
        AutoDiff (T aval) : val(aval)
        {
            for (size_t i = 0; i < D; i++) dval[i] = SCAL(0);
        }

        /// value with derivative 1 in direction diffindex; the int overload keeps AutoDiff(x, 0) unambiguous
        AutoDiff (SCAL aval, int diffindex) : AutoDiff(aval, size_t(diffindex)) { }

        AutoDiff (SCAL aval, size_t diffindex) : AutoDiff(aval)
        {
            dval[diffindex] = SCAL(1);
        }

        AutoDiff (SCAL aval, const SCAL * grad) : val(aval)
        {
            LoadGradient(grad);
        }

        /// evaluate an expression in one pass over the derivatives
        template <typename E, typename std::enable_if<is_node<E>, int>::type = 0>
        AutoDiff (const E & e)
        {
            static_assert(E::DIM == D, "AutoDiff dimensions must match");
            for (size_t i = 0; i < D; i++) dval[i] = e.DValue(i);
            val = e.Value();
        }

        template <typename E, typename std::enable_if<is_node<E>, int>::type = 0>
        AutoDiff & operator= (const E & e)
        {
            static_assert(E::DIM == D, "AutoDiff dimensions must match");
            // Nodes cached their values when built, and derivative i only reads index i
            // of its operands, so x = x*y + x is safe and the loop has no dependencies
            // between iterations. Saying so lets the compiler vectorize it even though
            // *this may appear on the right-hand side.
            const SCAL v = e.Value();
            ASC_BLA_AD_INDEPENDENT_ITERATIONS
            for (size_t i = 0; i < D; i++) dval[i] = e.DValue(i);
            val = v;
            return *this;
        }

        AutoDiff & operator= (SCAL aval)
        {
            val = aval;
            for (size_t i = 0; i < D; i++) dval[i] = SCAL(0);
            return *this;
        }

        template <typename T, typename std::enable_if<!ad_detail::is_expr<T> && std::is_convertible_v<T, SCAL>, int>::type = 0>
        AutoDiff & operator= (T aval)
        {
            val = aval;
            for (size_t i = 0; i < D; i++) dval[i] = SCAL(0);
            return *this;
        }

        void StoreGradient (SCAL * p) const { for (size_t i = 0; i < D; i++) p[i] = dval[i]; }
        void LoadGradient (const SCAL * p) { for (size_t i = 0; i < D; i++) dval[i] = p[i]; }

        SCAL Value() const { return val; }
        SCAL DValue (size_t i) const { return dval[i]; }
        SCAL * DValue () { return dval; }
        SCAL & Value() { return val; }
        SCAL & DValue (size_t i) { return dval[i]; }

        AutoDiff & operator+= (const AutoDiff & y)
        {
            val += y.val;
            for (size_t i = 0; i < D; i++) dval[i] += y.dval[i];
            return *this;
        }

        AutoDiff & operator-= (const AutoDiff & y)
        {
            val -= y.val;
            for (size_t i = 0; i < D; i++) dval[i] -= y.dval[i];
            return *this;
        }

        AutoDiff & operator*= (const AutoDiff & y)
        {
            const SCAL yv = y.val;
            for (size_t i = 0; i < D; i++) dval[i] = dval[i] * yv + val * y.dval[i];
            val *= yv;
            return *this;
        }

        AutoDiff & operator/= (const AutoDiff & y)
        {
            const SCAL inv = SCAL(1.0) / y.val;
            const SCAL q = val * inv;
            for (size_t i = 0; i < D; i++) dval[i] = (dval[i] - q * y.dval[i]) * inv;
            val = q;
            return *this;
        }

        template <typename E, typename std::enable_if<is_node<E>, int>::type = 0>
        AutoDiff & operator+= (const E & e)
        {
            const SCAL v = val + e.Value();
            ASC_BLA_AD_INDEPENDENT_ITERATIONS
            for (size_t i = 0; i < D; i++) dval[i] += e.DValue(i);
            val = v;
            return *this;
        }

        template <typename E, typename std::enable_if<is_node<E>, int>::type = 0>
        AutoDiff & operator-= (const E & e)
        {
            const SCAL v = val - e.Value();
            ASC_BLA_AD_INDEPENDENT_ITERATIONS
            for (size_t i = 0; i < D; i++) dval[i] -= e.DValue(i);
            val = v;
            return *this;
        }
        template <typename E, typename std::enable_if<is_node<E>, int>::type = 0>
        AutoDiff & operator*= (const E & e) { return *this *= AutoDiff(e); }
        template <typename E, typename std::enable_if<is_node<E>, int>::type = 0>
        AutoDiff & operator/= (const E & e) { return *this /= AutoDiff(e); }

        AutoDiff & operator*= (const SCAL & y)
        {
            for (size_t i = 0; i < D; i++) dval[i] *= y;
            val *= y;
            return *this;
        }

        AutoDiff & operator/= (const SCAL & y)
        {
            const SCAL iy = SCAL(1.0) / y;
            for (size_t i = 0; i < D; i++) dval[i] *= iy;
            val *= iy;
            return *this;
        }
    };

    template <size_t D, typename SCAL>
    inline std::ostream & operator<< (std::ostream & ost, const AutoDiff<D,SCAL> & x)
    {
        ost << x.Value() << ", D = ";
        for (size_t i = 0; i < D; i++)
            ost << x.DValue(i) << " ";
        return ost;
    }


    // ---------------------------------------------------------------------
    // Expression operators for AutoDiff<D, double> and its expressions
    // ---------------------------------------------------------------------
    // Defined in ad_detail, the namespace of the expression nodes and of
    // AutoDiff's base class, so argument-dependent lookup always finds them.
    namespace ad_detail {
        template <typename A, typename B>
        constexpr bool both_expr = is_expr<A> && is_expr<B>;


    template <typename A, typename B, typename std::enable_if<both_expr<A, B>, int>::type = 0>
    inline auto operator+ (A && a, B && b)
    { return Sum<Hold<A&&>, Hold<B&&>>(std::forward<A>(a), std::forward<B>(b)); }

    template <typename A, typename B, typename std::enable_if<both_expr<A, B>, int>::type = 0>
    inline auto operator- (A && a, B && b)
    { return Diff<Hold<A&&>, Hold<B&&>>(std::forward<A>(a), std::forward<B>(b)); }

    template <typename A, typename B, typename std::enable_if<both_expr<A, B>, int>::type = 0>
    inline auto operator* (A && a, B && b)
    { return Prod<Hold<A&&>, Hold<B&&>>(std::forward<A>(a), std::forward<B>(b)); }

    template <typename A, typename B, typename std::enable_if<both_expr<A, B>, int>::type = 0>
    inline auto operator/ (A && a, B && b)
    { return Quot<Hold<A&&>, Hold<B&&>>(std::forward<A>(a), std::forward<B>(b)); }

    template <typename A, typename std::enable_if<is_expr<A>, int>::type = 0>
    inline auto operator- (A && a)
    { return Affine<Hold<A&&>>(std::forward<A>(a), -1.0, 0.0); }

    template <typename S, typename A, typename std::enable_if<is_scalar<S> && is_expr<A>, int>::type = 0>
    inline auto operator* (S c, A && a)
    { return Scale<Hold<A&&>>(std::forward<A>(a), double(c)); }

    template <typename A, typename S, typename std::enable_if<is_scalar<S> && is_expr<A>, int>::type = 0>
    inline auto operator* (A && a, S c)
    { return Scale<Hold<A&&>>(std::forward<A>(a), double(c)); }

    template <typename A, typename S, typename std::enable_if<is_scalar<S> && is_expr<A>, int>::type = 0>
    inline auto operator/ (A && a, S c)
    { return Scale<Hold<A&&>>(std::forward<A>(a), 1.0 / double(c)); }

    template <typename S, typename A, typename std::enable_if<is_scalar<S> && is_expr<A>, int>::type = 0>
    inline auto operator/ (S c, A && a)
    { return Recip<Hold<A&&>>(std::forward<A>(a), double(c)); }

    template <typename S, typename A, typename std::enable_if<is_scalar<S> && is_expr<A>, int>::type = 0>
    inline auto operator+ (S c, A && a)
    { return Affine<Hold<A&&>>(std::forward<A>(a), 1.0, double(c)); }

    template <typename A, typename S, typename std::enable_if<is_scalar<S> && is_expr<A>, int>::type = 0>
    inline auto operator+ (A && a, S c)
    { return Affine<Hold<A&&>>(std::forward<A>(a), 1.0, double(c)); }

    template <typename S, typename A, typename std::enable_if<is_scalar<S> && is_expr<A>, int>::type = 0>
    inline auto operator- (S c, A && a)
    { return Affine<Hold<A&&>>(std::forward<A>(a), -1.0, double(c)); }

    template <typename A, typename S, typename std::enable_if<is_scalar<S> && is_expr<A>, int>::type = 0>
    inline auto operator- (A && a, S c)
    { return Affine<Hold<A&&>>(std::forward<A>(a), 1.0, -double(c)); }
    } // namespace ad_detail

    using ad_detail::operator+;
    using ad_detail::operator-;
    using ad_detail::operator*;
    using ad_detail::operator/;


    // ---------------------------------------------------------------------
    // Eager operators for nested types, AutoDiff<D, SCAL> with SCAL != double
    // ---------------------------------------------------------------------
    namespace ad_detail {
        template <typename SCAL>
        constexpr bool nested = !std::is_same_v<SCAL, double>;
    }

    template<size_t D, typename SCAL, typename std::enable_if<ad_detail::nested<SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator+ (const AutoDiff<D, SCAL> & x, const AutoDiff<D, SCAL> & y)
    {
        AutoDiff<D, SCAL> res;
        res.Value() = x.Value() + y.Value();
        for (size_t i = 0; i < D; i++) res.DValue(i) = x.DValue(i) + y.DValue(i);
        return res;
    }

    template<size_t D, typename SCAL, typename std::enable_if<ad_detail::nested<SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator- (const AutoDiff<D, SCAL> & x, const AutoDiff<D, SCAL> & y)
    {
        AutoDiff<D, SCAL> res;
        res.Value() = x.Value() - y.Value();
        for (size_t i = 0; i < D; i++) res.DValue(i) = x.DValue(i) - y.DValue(i);
        return res;
    }

    template<size_t D, typename SCAL, typename std::enable_if<ad_detail::nested<SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator- (const AutoDiff<D, SCAL> & x)
    {
        AutoDiff<D, SCAL> res;
        res.Value() = -x.Value();
        for (size_t i = 0; i < D; i++) res.DValue(i) = -x.DValue(i);
        return res;
    }

    template<size_t D, typename SCAL, typename std::enable_if<ad_detail::nested<SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator* (const AutoDiff<D, SCAL> & x, const AutoDiff<D, SCAL> & y)
    {
        AutoDiff<D, SCAL> res;
        const SCAL hx = x.Value(), hy = y.Value();
        res.Value() = hx * hy;
        for (size_t i = 0; i < D; i++) res.DValue(i) = hx * y.DValue(i) + hy * x.DValue(i);
        return res;
    }

    template<size_t D, typename SCAL, typename std::enable_if<ad_detail::nested<SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator/ (const AutoDiff<D, SCAL> & x, const AutoDiff<D, SCAL> & y)
    {
        AutoDiff<D, SCAL> res;
        const SCAL inv = 1.0 / y.Value();
        const SCAL q = x.Value() * inv;
        res.Value() = q;
        for (size_t i = 0; i < D; i++) res.DValue(i) = (x.DValue(i) - q * y.DValue(i)) * inv;
        return res;
    }

    template<size_t D, typename SCAL, typename SCAL2,
             typename std::enable_if<ad_detail::nested<SCAL> && std::is_convertible_v<SCAL2, SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator+ (SCAL2 x, const AutoDiff<D, SCAL> & y)
    {
        AutoDiff<D, SCAL> res = y;
        res.Value() = x + y.Value();
        return res;
    }

    template<size_t D, typename SCAL, typename SCAL2,
             typename std::enable_if<ad_detail::nested<SCAL> && std::is_convertible_v<SCAL2, SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator+ (const AutoDiff<D, SCAL> & y, SCAL2 x)
    {
        return x + y;
    }

    template<size_t D, typename SCAL, typename SCAL2,
             typename std::enable_if<ad_detail::nested<SCAL> && std::is_convertible_v<SCAL2, SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator- (const AutoDiff<D, SCAL> & x, SCAL2 y)
    {
        AutoDiff<D, SCAL> res = x;
        res.Value() = x.Value() - y;
        return res;
    }

    template<size_t D, typename SCAL, typename SCAL2,
             typename std::enable_if<ad_detail::nested<SCAL> && std::is_convertible_v<SCAL2, SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator- (SCAL2 x, const AutoDiff<D, SCAL> & y)
    {
        AutoDiff<D, SCAL> res;
        res.Value() = x - y.Value();
        for (size_t i = 0; i < D; i++) res.DValue(i) = -y.DValue(i);
        return res;
    }

    template<size_t D, typename SCAL, typename SCAL2,
             typename std::enable_if<ad_detail::nested<SCAL> && std::is_convertible_v<SCAL2, SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator* (const SCAL2 x, const AutoDiff<D, SCAL> & y)
    {
        AutoDiff<D, SCAL> res;
        res.Value() = x * y.Value();
        for (size_t i = 0; i < D; i++) res.DValue(i) = x * y.DValue(i);
        return res;
    }

    template<size_t D, typename SCAL, typename SCAL2,
             typename std::enable_if<ad_detail::nested<SCAL> && std::is_convertible_v<SCAL2, SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator* (const AutoDiff<D, SCAL> & y, SCAL2 x)
    {
        return x * y;
    }

    template<size_t D, typename SCAL, typename SCAL2,
             typename std::enable_if<ad_detail::nested<SCAL> && std::is_convertible_v<SCAL2, SCAL>, int>::type = 0>
    inline AutoDiff<D, SCAL> operator/ (const AutoDiff<D, SCAL> & x, SCAL2 y)
    {
        return (1.0 / y) * x;
    }


    // ---------------------------------------------------------------------
    // Functions. Expression arguments are evaluated first.
    // ---------------------------------------------------------------------
    template<size_t D, typename SCAL>
    inline AutoDiff<D, SCAL> Inv (const AutoDiff<D, SCAL> & x)
    {
        AutoDiff<D, SCAL> res;
        const SCAL inv = 1.0 / x.Value();
        res.Value() = inv;
        const SCAL f = -inv * inv;
        for (size_t i = 0; i < D; i++) res.DValue(i) = f * x.DValue(i);
        return res;
    }

    using std::sqrt;
    template<size_t D, typename SCAL>
    inline AutoDiff<D, SCAL> sqrt (const AutoDiff<D, SCAL> & x)
    {
        using std::sqrt;
        AutoDiff<D, SCAL> res;
        res.Value() = sqrt(x.Value());
        const SCAL f = 0.5 / res.Value();
        for (size_t j = 0; j < D; j++) res.DValue(j) = f * x.DValue(j);
        return res;
    }

    template <size_t D, typename SCAL>
    inline AutoDiff<D, SCAL> exp (const AutoDiff<D, SCAL> & x)
    {
        using std::exp;
        AutoDiff<D, SCAL> res;
        res.Value() = exp(x.Value());
        for (size_t k = 0; k < D; k++) res.DValue(k) = x.DValue(k) * res.Value();
        return res;
    }

    template <size_t D, typename SCAL>
    inline AutoDiff<D, SCAL> log (const AutoDiff<D, SCAL> & x)
    {
        using std::log;
        AutoDiff<D, SCAL> res;
        res.Value() = log(x.Value());
        const SCAL xinv = 1.0 / x.Value();
        for (size_t k = 0; k < D; k++) res.DValue(k) = x.DValue(k) * xinv;
        return res;
    }

    using std::pow;
    template <size_t D, typename SCAL>
    inline AutoDiff<D,SCAL> pow (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y)
    {
        return exp(AutoDiff<D,SCAL>(log(x) * y));
    }

    template <size_t D, typename SCAL>
    inline AutoDiff<D, SCAL> pow (const AutoDiff<D, SCAL> & x, SCAL y)
    {
        // Direct power rule, so negative and zero bases work for real exponents
        // where the result is defined (exp(y*log(x)) would give NaN there).
        using std::pow;
        AutoDiff<D, SCAL> res;
        res.Value() = pow(x.Value(), y);
        const SCAL factor = (y == SCAL(0)) ? SCAL(0) : y * pow(x.Value(), y - SCAL(1));
        for (size_t k = 0; k < D; k++) res.DValue(k) = factor * x.DValue(k);
        return res;
    }

    template <size_t D, typename SCAL>
    inline AutoDiff<D, SCAL> sin (const AutoDiff<D, SCAL> & x)
    {
        using std::sin; using std::cos;
        AutoDiff<D, SCAL> res;
        const SCAL c = cos(x.Value());
        res.Value() = sin(x.Value());
        for (size_t k = 0; k < D; k++) res.DValue(k) = x.DValue(k) * c;
        return res;
    }

    template <size_t D, typename SCAL>
    inline AutoDiff<D, SCAL> cos (const AutoDiff<D, SCAL> & x)
    {
        using std::sin; using std::cos;
        AutoDiff<D, SCAL> res;
        const SCAL s = sin(x.Value());
        res.Value() = cos(x.Value());
        for (size_t k = 0; k < D; k++) res.DValue(k) = -s * x.DValue(k);
        return res;
    }

    template <size_t D, typename SCAL>
    inline AutoDiff<D, SCAL> tan (const AutoDiff<D, SCAL> & x)
    {
        using std::tan;
        AutoDiff<D, SCAL> res;
        const SCAL t = tan(x.Value());
        res.Value() = t;
        const SCAL f = 1.0 + t * t;
        for (size_t k = 0; k < D; k++) res.DValue(k) = f * x.DValue(k);
        return res;
    }

#define ASC_BLA_AD_EXPR_FUNCTION(NAME) \
    template <typename E, typename std::enable_if<ad_detail::is_node<E>, int>::type = 0> \
    inline auto NAME (const E & e) { return NAME(AutoDiff<E::DIM, double>(e)); }
    ASC_BLA_AD_EXPR_FUNCTION(Inv)
    ASC_BLA_AD_EXPR_FUNCTION(sqrt)
    ASC_BLA_AD_EXPR_FUNCTION(exp)
    ASC_BLA_AD_EXPR_FUNCTION(log)
    ASC_BLA_AD_EXPR_FUNCTION(sin)
    ASC_BLA_AD_EXPR_FUNCTION(cos)
    ASC_BLA_AD_EXPR_FUNCTION(tan)
#undef ASC_BLA_AD_EXPR_FUNCTION

    template <typename E, typename std::enable_if<ad_detail::is_node<E>, int>::type = 0>
    inline auto pow (const E & e, double y) { return pow(AutoDiff<E::DIM, double>(e), y); }
}

#endif
