#include <bairstow/ThreadPool.h> // for ThreadPool

#include <cstddef> // for size_t

#include <bairstow/rootfinding.hpp> // for Vec2, delta, Options, horner_eval
#include <cmath>                    // for abs, acos, cos, pow
#include <functional>               // for __base
#include <future>                   // for future
#include <py2cpp/range.hpp>         // for range
#include <thread>                   // for thread
#include <type_traits>              // for move
#include <utility>                  // for pair
#include <vector>                   // for vector, vector<>::reference, __v...

#include "bairstow/vector2.hpp" // for operator-, Vector2

// using Vec2 = numeric::Vector2<double>;
// using Mat2 = numeric::Matrix2<Vec2>;

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] vr
 * @return Vec2
 */
auto horner(std::vector<double> &pb, size_t n, const Vec2 &vr) -> Vec2 {
  const auto &r = vr.x();
  const auto &t = vr.y();
  pb[1] -= pb[0] * r;
  for (auto i : py::range(2, n)) {
    pb[i] -= pb[i - 1] * r + pb[i - 2] * t;
  }
  pb[n] -= pb[n - 2] * t;
  return Vec2{pb[n - 1], pb[n]};
}

/**
 * @brief
 *
 * @param[in] pa
 * @return std::vector<Vec2>
 */
auto initial_guess(const std::vector<double> &pa) -> std::vector<Vec2> {
  static const auto PI = std::acos(-1.);

  auto N = int(pa.size()) - 1;
  const auto c = -pa[1] / (N * pa[0]);
  auto pb = pa;
  const auto Pc = horner_eval(pb, N, c); // TODO
  const auto re = std::pow(std::abs(Pc), 1. / N);
  N /= 2;
  N *= 2; // make even
  const auto k = PI / N;
  const auto m = c * c + re * re;
  auto vr0s = std::vector<Vec2>{};
  for (auto i = 1; i < N; i += 2) {
    const auto temp = re * std::cos(k * i);
    auto r0 = -2 * (c + temp);
    auto t0 = m + 2 * c * temp;
    vr0s.emplace_back(Vec2{std::move(r0), std::move(t0)});
  }
  return vr0s;
}

/**
 * @brief Multi-threading Bairstow's method (even degree only)
 *
 * @param[in] pa polynomial
 * @param[in,out] vrs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return std::pair<unsigned int, bool>
 */
auto pbairstow_even(const std::vector<double> &pa, std::vector<Vec2> &vrs,
                    const Options &options = Options())
    -> std::pair<unsigned int, bool> {
  const auto N = pa.size() - 1; // degree, assume even
  const auto M = vrs.size();
  auto found = false;
  auto converged = std::vector<bool>(M, false);
  auto niter = 1U;
  ThreadPool pool(std::thread::hardware_concurrency());

  for (; niter != options.max_iter; ++niter) {
    auto tol = 0.0;
    std::vector<std::future<double>> results;

    for (auto i : py::range(M)) {
      if (converged[i]) {
        continue;
      }
      results.emplace_back(pool.enqueue([&, i]() {
        auto pb = pa;
        // auto n = pa.size() - 1;
        const auto &vri = vrs[i];
        const auto vA = horner(pb, N, vri);
        const auto tol_i = std::max(std::abs(vA.x()), std::abs(vA.y()));
        if (tol_i < 1e-15) {
          converged[i] = true;
          return tol_i;
        }
        auto vA1 = horner(pb, N - 2, vri);
        for (auto j : py::range(M)) { // exclude i
          if (j == i) {
            continue;
          }
          const auto vrj = vrs[j]; // make a copy, don't reference!
          vA1 -= delta(vA, vrj, vri - vrj);
        }
        vrs[i] -= delta(vA, vri, std::move(vA1)); // Gauss-Seidel fashion
        return tol_i;
      }));
    }
    for (auto &&result : results) {
      auto &&res = result.get();
      if (tol < res) {
        tol = res;
      }
    }
    if (tol < options.tol) {
      found = true;
      break;
    }
  }
  return {niter, found};
}

// auto find_rootq(const Vec2& r) {
//     auto hb = b / 2.;
//     auto d = hb * hb - c;
//     if (d < 0.) {
//         auto x1 = -hb + (sqrt(-d) if (hb < 0. else -sqrt(-d))*1j;
//     }
//     else {
//         auto x1 = -hb + (sqrt(d) if (hb < 0. else -sqrt(d));
//     }
//     auto x2 = c / x1;
//     return x1, x2;
// }
