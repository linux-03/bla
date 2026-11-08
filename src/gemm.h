#ifndef FILE_GEMM_H
#define FILE_GEMM_H

// Matrix product C = alpha * A * B + beta * C for double precision, without BLAS.
//
// The method follows Goto and van de Geijn, "Anatomy of High-Performance Matrix
// Multiplication" (2008), as structured in BLIS:
//   * B is cut into KC x NC blocks and A into MC x KC blocks, sized for the caches.
//   * Each block is packed into a contiguous buffer in the order the kernel reads it,
//     so any layout (row-major, column-major, transposed, sub-view) is read sequentially.
//   * A micro-kernel keeps an MR x NR block of C in SIMD registers and updates it with
//     one fused multiply-add per register and step.
//   * Large products split the MC blocks across threads.

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <vector>

#include "simd.h"

namespace ASC_bla {

  namespace gemm_detail {

    // Register block: MR rows of C, each held as one SIMD<double, NR>.
    // 6 x 8 uses 24 of 32 NEON registers (or 12 of 16 AVX registers) for C.
    constexpr size_t MR = 6;
    constexpr size_t NR = 8;
    // Cache blocks. KC x NR of B (16 KB) stays in L1, MC x KC of A (240 KB) in L2.
    constexpr size_t KC = 256;
    constexpr size_t MC = 120;
    constexpr size_t NC = 4096;

    using Row = ASC_HPC::SIMD<double, NR>;

    inline std::atomic<int> & ThreadSetting()
    {
#ifdef __EMSCRIPTEN__
      static std::atomic<int> threads{1};  // no threads in a plain WebAssembly build
#else
      static std::atomic<int> threads{int(std::max(1u, std::thread::hardware_concurrency()))};
#endif
      return threads;
    }

    // Pack an mc x kc block of A into panels of MR rows: panel, then k, then row.
    inline void PackA(size_t mc, size_t kc, const double * a, ptrdiff_t rsa, ptrdiff_t csa, double * buf)
    {
      for (size_t i0 = 0; i0 < mc; i0 += MR)
        {
          const size_t rows = std::min(MR, mc - i0);
          for (size_t p = 0; p < kc; p++)
            {
              const double * col = a + ptrdiff_t(i0) * rsa + ptrdiff_t(p) * csa;
              size_t i = 0;
              for ( ; i < rows; i++) buf[i] = col[ptrdiff_t(i) * rsa];
              for ( ; i < MR; i++) buf[i] = 0.0;
              buf += MR;
            }
        }
    }

    // Pack a kc x nc block of B into panels of NR columns: panel, then k, then column.
    inline void PackB(size_t kc, size_t nc, const double * b, ptrdiff_t rsb, ptrdiff_t csb, double * buf)
    {
      for (size_t j0 = 0; j0 < nc; j0 += NR)
        {
          const size_t cols = std::min(NR, nc - j0);
          for (size_t p = 0; p < kc; p++)
            {
              const double * row = b + ptrdiff_t(p) * rsb + ptrdiff_t(j0) * csb;
              size_t j = 0;
              if (csb == 1 && cols == NR)
                for ( ; j < NR; j++) buf[j] = row[j];
              else
                {
                  for ( ; j < cols; j++) buf[j] = row[ptrdiff_t(j) * csb];
                  for ( ; j < NR; j++) buf[j] = 0.0;
                }
              buf += NR;
            }
        }
    }

    // C[0:rows, 0:cols] += alpha * (packed A panel) * (packed B panel).
    inline void MicroKernel(size_t kc, const double * ap, const double * bp, double alpha,
                            double * c, ptrdiff_t rsc, ptrdiff_t csc, size_t rows, size_t cols)
    {
      Row acc0(0.0), acc1(0.0), acc2(0.0), acc3(0.0), acc4(0.0), acc5(0.0);
      for (size_t p = 0; p < kc; p++)
        {
          Row b(const_cast<double*>(bp));
          acc0 = FMA(Row(ap[0]), b, acc0);
          acc1 = FMA(Row(ap[1]), b, acc1);
          acc2 = FMA(Row(ap[2]), b, acc2);
          acc3 = FMA(Row(ap[3]), b, acc3);
          acc4 = FMA(Row(ap[4]), b, acc4);
          acc5 = FMA(Row(ap[5]), b, acc5);
          ap += MR;
          bp += NR;
        }

      const Row scale(alpha);
      Row * acc[MR] = { &acc0, &acc1, &acc2, &acc3, &acc4, &acc5 };
      if (rows == MR && cols == NR && csc == 1)
        {
          for (size_t i = 0; i < MR; i++)
            {
              double * ci = c + ptrdiff_t(i) * rsc;
              FMA(scale, *acc[i], Row(ci)).Store(ci);
            }
          return;
        }

      double tile[MR * NR];
      for (size_t i = 0; i < MR; i++) acc[i]->Store(tile + i * NR);
      for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
          c[ptrdiff_t(i) * rsc + ptrdiff_t(j) * csc] += alpha * tile[i * NR + j];
    }

    // Like MicroKernel, but reads A and B in place: B rows must be contiguous and
    // hold at least NR columns. Used for small products, where packing does not pay off.
    inline void DirectKernel(size_t kc, const double * a, ptrdiff_t rsa, ptrdiff_t csa,
                             const double * b, ptrdiff_t rsb, double alpha,
                             double * c, ptrdiff_t rsc, ptrdiff_t csc, size_t rows)
    {
      // Rows past the end repeat the last row; their results are discarded below.
      // Named accumulators (not an array) let the compiler keep them in registers.
      const double * a0 = a;
      const double * a1 = a + ptrdiff_t(std::min<size_t>(1, rows - 1)) * rsa;
      const double * a2 = a + ptrdiff_t(std::min<size_t>(2, rows - 1)) * rsa;
      const double * a3 = a + ptrdiff_t(std::min<size_t>(3, rows - 1)) * rsa;
      const double * a4 = a + ptrdiff_t(std::min<size_t>(4, rows - 1)) * rsa;
      const double * a5 = a + ptrdiff_t(std::min<size_t>(5, rows - 1)) * rsa;
      Row acc0(0.0), acc1(0.0), acc2(0.0), acc3(0.0), acc4(0.0), acc5(0.0);
      for (size_t p = 0; p < kc; p++)
        {
          Row bp(const_cast<double*>(b + ptrdiff_t(p) * rsb));
          const ptrdiff_t q = ptrdiff_t(p) * csa;
          acc0 = FMA(Row(a0[q]), bp, acc0);
          acc1 = FMA(Row(a1[q]), bp, acc1);
          acc2 = FMA(Row(a2[q]), bp, acc2);
          acc3 = FMA(Row(a3[q]), bp, acc3);
          acc4 = FMA(Row(a4[q]), bp, acc4);
          acc5 = FMA(Row(a5[q]), bp, acc5);
        }

      double tile[MR * NR];
      acc0.Store(tile);          acc1.Store(tile + NR);     acc2.Store(tile + 2 * NR);
      acc3.Store(tile + 3 * NR); acc4.Store(tile + 4 * NR); acc5.Store(tile + 5 * NR);
      for (size_t i = 0; i < rows; i++)
        {
          double * ci = c + ptrdiff_t(i) * rsc;
          if (csc == 1)
            FMA(Row(alpha), Row(tile + i * NR), Row(ci)).Store(ci);
          else
            for (size_t j = 0; j < NR; j++) ci[ptrdiff_t(j) * csc] += alpha * tile[i * NR + j];
        }
    }

    // All register blocks of one packed MC x KC block of A against the packed B block.
    inline void MacroKernel(size_t mc, size_t nc, size_t kc, const double * apack, const double * bpack,
                            double alpha, double * c, ptrdiff_t rsc, ptrdiff_t csc)
    {
      for (size_t j0 = 0; j0 < nc; j0 += NR)
        {
          const size_t cols = std::min(NR, nc - j0);
          const double * bp = bpack + j0 * kc;
          for (size_t i0 = 0; i0 < mc; i0 += MR)
            MicroKernel(kc, apack + i0 * kc, bp, alpha,
                        c + ptrdiff_t(i0) * rsc + ptrdiff_t(j0) * csc, rsc, csc,
                        std::min(MR, mc - i0), cols);
        }
    }

    inline size_t RoundUp(size_t x, size_t to) { return (x + to - 1) / to * to; }

    // Reusable packing buffer for the calling thread, so small products do not allocate.
    inline double * ScratchBuffer(int which, size_t size)
    {
      thread_local std::vector<double> buffers[2];
      if (buffers[which].size() < size) buffers[which].resize(size);
      return buffers[which].data();
    }

  } // namespace gemm_detail


  // Number of threads used for large products. Defaults to all hardware threads.
  inline void SetGemmThreads(int threads) { gemm_detail::ThreadSetting() = std::max(1, threads); }
  inline int GetGemmThreads() { return gemm_detail::ThreadSetting(); }


  // C = alpha * A * B + beta * C, where A is m x k, B is k x n and C is m x n.
  // Element (i, j) of X is at x[i * rsx + j * csx]: any row-major, column-major,
  // transposed or sub-matrix layout. C must not overlap A or B.
  inline void Gemm(size_t m, size_t n, size_t k, double alpha,
                   const double * a, ptrdiff_t rsa, ptrdiff_t csa,
                   const double * b, ptrdiff_t rsb, ptrdiff_t csb,
                   double beta, double * c, ptrdiff_t rsc, ptrdiff_t csc)
  {
    using namespace gemm_detail;
    if (m == 0 || n == 0) return;

    // beta == 0 overwrites C, so NaN or garbage in C does not propagate.
    if (beta != 1.0)
      for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
          {
            double & cij = c[ptrdiff_t(i) * rsc + ptrdiff_t(j) * csc];
            cij = (beta == 0.0) ? 0.0 : beta * cij;
          }
    if (k == 0 || alpha == 0.0) return;

    // Tiny products: packing would cost more than it saves.
    if (n < NR || m * n * k <= 64)
      {
        for (size_t i = 0; i < m; i++)
          for (size_t j = 0; j < n; j++)
            {
              double s = 0;
              for (size_t p = 0; p < k; p++)
                s += a[ptrdiff_t(i) * rsa + ptrdiff_t(p) * csa] * b[ptrdiff_t(p) * rsb + ptrdiff_t(j) * csb];
              c[ptrdiff_t(i) * rsc + ptrdiff_t(j) * csc] += alpha * s;
            }
        return;
      }

    // Small products with contiguous rows of B: read B in place instead of packing.
    // Everything fits in L1 here, so packing would only add work.
    if (csb == 1 && k <= KC && m <= MC && n <= 64)
      {
        for (size_t j0 = 0; j0 + NR <= n; j0 += NR)
          for (size_t i0 = 0; i0 < m; i0 += MR)
            DirectKernel(k, a + ptrdiff_t(i0) * rsa, rsa, csa, b + j0, rsb, alpha,
                         c + ptrdiff_t(i0) * rsc + ptrdiff_t(j0) * csc, rsc, csc, std::min(MR, m - i0));
        const size_t done = n / NR * NR;
        for (size_t i = 0; i < m; i++)
          for (size_t j = done; j < n; j++)
            {
              double s = 0;
              for (size_t p = 0; p < k; p++)
                s += a[ptrdiff_t(i) * rsa + ptrdiff_t(p) * csa] * b[ptrdiff_t(p) * rsb + ptrdiff_t(j)];
              c[ptrdiff_t(i) * rsc + ptrdiff_t(j) * csc] += alpha * s;
            }
        return;
      }

    const size_t mblocks = (m + MC - 1) / MC;
    // Threads pay off only when there is enough work per thread.
    size_t threads = 1;
    if (double(m) * double(n) * double(k) >= 192.0 * 192.0 * 192.0)
      threads = std::min<size_t>(size_t(GetGemmThreads()), mblocks);

    const size_t asize = RoundUp(std::min(m, MC), MR) * KC;
    double * bpack = ScratchBuffer(0, KC * RoundUp(std::min(n, NC), NR));
    std::vector<double *> apacks(threads);
    std::vector<std::vector<double>> extra(threads > 1 ? threads - 1 : 0, std::vector<double>(asize));
    apacks[0] = ScratchBuffer(1, asize);
    for (size_t t = 1; t < threads; t++) apacks[t] = extra[t - 1].data();

    for (size_t jc = 0; jc < n; jc += NC)
      {
        const size_t nc = std::min(NC, n - jc);
        for (size_t pc = 0; pc < k; pc += KC)
          {
            const size_t kc = std::min(KC, k - pc);
            PackB(kc, nc, b + ptrdiff_t(pc) * rsb + ptrdiff_t(jc) * csb, rsb, csb, bpack);

            auto work = [&](size_t worker)
            {
              double * apack = apacks[worker];
              for (size_t block = worker; block < mblocks; block += threads)
                {
                  const size_t ic = block * MC;
                  const size_t mc = std::min(MC, m - ic);
                  PackA(mc, kc, a + ptrdiff_t(ic) * rsa + ptrdiff_t(pc) * csa, rsa, csa, apack);
                  MacroKernel(mc, nc, kc, apack, bpack, alpha,
                              c + ptrdiff_t(ic) * rsc + ptrdiff_t(jc) * csc, rsc, csc);
                }
            };

            if (threads == 1)
              work(0);
            else
              {
                std::vector<std::thread> pool;
                for (size_t t = 1; t < threads; t++) pool.emplace_back(work, t);
                work(0);
                for (auto & t : pool) t.join();
              }
          }
      }
  }

} // namespace ASC_bla

#endif
