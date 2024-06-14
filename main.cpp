#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <random>
#include <thread>
#include <unordered_set>
#include <vector>

#include "json.hpp"

using namespace std;
using namespace std::chrono;
using namespace nlohmann;

/**
 * A UBQP problem instance.
 */
struct ubqp {
  /** The number of variables in a solution vector. */
  size_t n;

  /** The matrix, stored in row-major order. */
  unique_ptr<long[]> Q;

  /**
   * Constructs an UBQP object.
   */
  ubqp(size_t n) : n(n), Q(new long[n * n]()) {}

  /**
   * Obtains the pointer to the beginning of a row of the matrix.
   */
  long* operator[](size_t i) const { return &Q[n * i]; }

  /**
   * Constructs a UBQP problem instance from a file.
   */
  static ubqp load(string label) {
    size_t n, number;
    if (sscanf(label.c_str(), "bqp%zu.%zu", &n, &number) == 2) return from_bqp(n, number);
    if (sscanf(label.c_str(), "G%zu", &number) == 1) return from_maxcut(number);
    assert(false);
    exit(1);
  }

  /**
   * Constructs a UBQP problem instance from an OR-Library file.
   */
  static ubqp from_bqp(size_t n, size_t number) {
    ifstream file("instances/bqp" + to_string(n) + string(".txt"), ifstream::in);
    assert(file);
    size_t count;
    file >> count;
    for (size_t k = 1; k <= count; k++) {
      size_t nonzeros, i, j;
      long q;
      if (k == number) {
        file >> n >> nonzeros;
        ubqp Q(n);
        for (size_t l = 0; l < nonzeros; l++) {
          file >> i >> j >> q;
          i -= 1;
          j -= 1;
          Q[i][j] = Q[j][i] = i == j ? -q : 2 * -q;
        }
        return Q;
      } else {
        file >> n >> nonzeros;
        for (size_t l = 0; l < nonzeros; l++) file >> i >> j >> q;
      }
    }
    assert(false);
    exit(1);
  }

  /**
   * Constructs a UBQP problem instance from a Max-Cut problem file.
   */
  static ubqp from_maxcut(size_t number) {
    ifstream file("instances/G" + to_string(number), ifstream::in);
    assert(file);
    size_t n, nonzeros;
    file >> n >> nonzeros;
    ubqp Q(n);
    for (size_t l = 0; l < nonzeros; l++) {
      size_t i, j;
      long q;
      file >> i >> j >> q;
      i -= 1;
      j -= 1;
      Q[i][i] -= q;
      Q[j][j] -= q;
      Q[i][j] += 2 * q;
      Q[j][i] += 2 * q;
    }
    return Q;
  }
};

/**
 * Represents a solution to a UBQP problem instance.
 */
struct solution {
  /** The objective function value. */
  long fx;
  /** The solution vector. */
  unique_ptr<long[]> x;
  /** The reevaluation vector. */
  unique_ptr<long[]> dx;

  /** Constructs a solution. */
  solution(const ubqp& Q) : fx(0), x(new long[Q.n]()), dx(new long[Q.n]()) {}

  /** Returns whether the solution is better than the other. */
  bool operator<(const solution& other) const { return fx < other.fx; }

  /** Computes a hash of the solution vector. */
  size_t hash(const ubqp& Q) const {
    size_t seed = 0;
    for (size_t i = 0; i < Q.n; i++) seed = seed * 31 + 2 * (x[i] != 0) + 3;
    return seed;
  }
};

/**
 * Copies a solution into another.
 */
void copy_solution(const ubqp& Q, const solution& x, solution& y) {
  copy_n(&x.x[0], Q.n, &y.x[0]);
  copy_n(&x.dx[0], Q.n, &y.dx[0]);
  y.fx = x.fx;
}

/**
 * Represents an evaluation algorithm to use when evaluating UBQP solutions.
 */
enum struct evaluation { basic, f_inc, rfflip, rfflip_rv };

/**
 * Checks whether the solution, when an incumbent, has the correct objective function value.
 */
bool check_incumbent_evaluation(const ubqp& Q, const solution& x) {
  long fx_tmp = 0;
  for (size_t i = 0; i < Q.n; i++)
    for (size_t j = 0; j <= i; j++) fx_tmp += x.x[i] * x.x[j] * Q[i][j];
  return x.fx == fx_tmp;
}

/**
 * Checks whether the solution, when a neighbor, has the correct objective function value.
 *
 * @param Q the problem instance.
 * @param x the incumbent solution.
 * @param y the neighbor solution to check.
 * @param r the number of flip moves between x and y.
 * @param N the permutation vector whose first r indices have been changed between x and y.
 */
bool check_evaluation(const ubqp& Q, const solution& x, solution& y, size_t r,
                      const unique_ptr<size_t[]>& N) {
  for (size_t k = r; k < Q.n; k++) y.x[N[k]] = x.x[N[k]];
  long fy_tmp = 0;
  for (size_t i = 0; i < Q.n; i++)
    for (size_t j = 0; j <= i; j++) fy_tmp += y.x[i] * y.x[j] * Q[i][j];
  return y.fx == fy_tmp;
}

/**
 * Checks whether the solution has a correct reevaluation vector.
 */
bool check_rv(const ubqp& Q, const solution& x) {
  for (size_t i = 0; i < Q.n; i++) {
    long dxi_tmp = 0;
    for (size_t j = 0; j < Q.n; j++)
      if (j != i) dxi_tmp += x.x[j] * Q[i][j];
    if (x.dx[i] != dxi_tmp) return false;
  }
  return true;
}

/**
 * Builds the reevaluation vector of a solution from scratch.
 */
void build_rv(const ubqp& Q, solution& x) {
  for (size_t i = 0; i < Q.n; i++) {
    x.dx[i] = 0;
    for (size_t j = 0; j < Q.n; j++)
      if (j != i) x.dx[i] += x.x[j] * Q[i][j];
  }
}

/**
 * Updates an incumbent's reevaluation vector by using a neighbor solution.
 */
void update_rv(const ubqp& Q, solution& x, const solution& y, size_t r,
               const unique_ptr<size_t[]>& N) {
  for (size_t i = 0; i < Q.n; i++)
    for (size_t m = 0; m < r; m++) {
      size_t j = N[m];
      if (i != j) x.dx[i] += Q[i][j] * (y.x[j] - x.x[j]);
    }
}

/**
 * Evaluates a one-flip binary move.
 *
 * @param eval the evaluation algorithm to use.
 * @param Q the problem instance.
 * @param x the solution.
 * @param i the index of the component to flip.
 * @return the difference in objective function value after the flip.
 */
long oneflip_evaluate(const ubqp& Q, const solution& x, size_t i) {
  return (1 - 2 * x.x[i]) * (x.dx[i] + Q[i][i]);
}

/**
 * Performs a one-flip binary move.
 *
 * @param eval the evaluation algorithm to use.
 * @param Q the problem instance.
 * @param cur the solution.
 * @param i the index of the component to flip.
 */
void oneflip_move(const ubqp& Q, solution& cur, size_t i) {
  for (size_t j = 0; j < i; j++) cur.dx[j] += Q[j][i] * (1 - 2 * cur.x[i]);
  for (size_t j = i + 1; j < Q.n; j++) cur.dx[j] += Q[j][i] * (1 - 2 * cur.x[i]);
  cur.fx = cur.fx + oneflip_evaluate(Q, cur, i);
  cur.x[i] = 1 - cur.x[i];
}

bool oneflip_move_greedy(const ubqp& Q, solution& cur, const unique_ptr<size_t[]>& N) {
  (void)N;
  size_t i = Q.n;
  for (size_t j = 0; j < Q.n; j++) {
    long obj = oneflip_evaluate(Q, cur, j);
    if (obj < 0 && (i == Q.n || obj < oneflip_evaluate(Q, cur, i))) i = j;
  }
  if (i == Q.n) return false;
  oneflip_move(Q, cur, i);
  return true;
}

void oneflip_move_greedy_perturb(const ubqp& Q, solution& cur, const unique_ptr<size_t[]>& N,
                                 size_t r) {
  (void)N;
  for (size_t k = 0; k < r; k++) {
    size_t i = Q.n;
    for (size_t j = 0; j < Q.n; j++) {
      long obj = oneflip_evaluate(Q, cur, j);
      if (obj >= 0 && (i == Q.n || obj < oneflip_evaluate(Q, cur, i))) i = j;
    }
    if (i == Q.n) return;
    oneflip_move(Q, cur, i);
  }
}

bool rflip_move_greedy(const ubqp& Q, solution& cur, const unique_ptr<size_t[]>& N) {
  size_t r =
      partition(&N[0], &N[Q.n], [&](size_t i) { return oneflip_evaluate(Q, cur, i) < 0; }) - &N[0];
  if (!r) return false;
  // for (size_t k = 0; k < r; k++) cout << oneflip_evaluate(Q, cur, N[k]) << " ";
  // cout << endl;
  for (size_t k = 0; k < r; k++) oneflip_move(Q, cur, N[k]);
  return true;
}

void rflip_move_greedy_perturb(const ubqp& Q, solution& cur, const unique_ptr<size_t[]>& N,
                               size_t r) {
  iota(&N[0], &N[Q.n], size_t(0));
  sort(&N[0], &N[Q.n], [&](size_t i, size_t j) {
    return oneflip_evaluate(Q, cur, i) < oneflip_evaluate(Q, cur, j);
  });
  for (size_t k = 0; k < r; k++) oneflip_move(Q, cur, N[k]);
}

/**
 * Performs a local improvement procedure.
 *
 * @param Q the problem instance.
 * @param x the incumbent solution to perform the search on.
 * @param y buffer for neighbor solutions.
 * @param F the maximum value within the components of x and y.
 * @param r the number of flip moves between x and y.
 * @param l_max the maximum amount of iterations without improvements.
 * @param L the vector for storing tabu moves.
 * @param N buffer for permutation vectors.
 * @param K the tabu tenure.
 * @param rng the random number generator.
 */
void ts(const ubqp& Q, solution& cur, solution& tmp, long F, size_t l_max, unique_ptr<size_t[]>& L,
        unique_ptr<size_t[]>& N, size_t K, mt19937& rng) {
  (void)F, (void)l_max, (void)L, (void)K;
  for (size_t c = 0; c < Q.n - 1; c += 1 + (rng() % (Q.n / 10))) {
    copy_solution(Q, cur, tmp);
    // oneflip_move_greedy_perturb(Q, tmp, N, c);
    rflip_move_greedy_perturb(Q, tmp, N, c);
    while (oneflip_move_greedy(Q, tmp, N)) continue;
    // while (rflip_move_greedy(Q, tmp, N)) continue;
    if (tmp < cur) {
      swap(tmp, cur);
      c = 0;
    }
  }
}

/**
 * Evaluates a neighbor solution.
 *
 * @param eval the evaluation algorithm to use.
 * @param Q the problem instance.
 * @param x the incumbent solution.
 * @param y the neighbor solution to evaluate.
 * @param r the number of flip moves between x and y.
 * @param N the permutation vector whose first r indices have been changed between x and y.
 */
template <evaluation eval>
void evaluate(const ubqp& Q, const solution& x, solution& y, size_t r,
              const unique_ptr<size_t[]>& N);

template <>
void evaluate<evaluation::basic>(const ubqp& Q, const solution& x, solution& y, size_t r,
                                 const unique_ptr<size_t[]>& N) {
  for (size_t k = r; k < Q.n; k++) y.x[N[k]] = x.x[N[k]];
  y.fx = 0;
  for (size_t i = 0; i < Q.n; i++)
    if (y.x[i])
      for (size_t l = 0; l <= i; l++) y.fx += y.x[i] * y.x[l] * Q[i][l];
}

template <>
void evaluate<evaluation::f_inc>(const ubqp& Q, const solution& x, solution& y, size_t r,
                                 const unique_ptr<size_t[]>& N) {
  for (size_t i = 0; i < Q.n; i++) y.dx[i] = x.dx[i];
  y.fx = x.fx;
  for (size_t m = 0; m < r; m++) {
    size_t i = N[m];
    y.fx += (y.x[i] - x.x[i]) * (y.dx[i] + Q[i][i] * (y.x[i] + x.x[i]));
    if (m < r - 1)
      for (size_t j = 0; j < Q.n; j++)
        if (j != i) y.dx[j] += Q[i][j] * (y.x[i] - x.x[i]);
  }
}

template <>
void evaluate<evaluation::rfflip>(const ubqp& Q, const solution& x, solution& y, size_t r,
                                  const unique_ptr<size_t[]>& N) {
  y.fx = x.fx;
  for (size_t m = 0; m < r; m++) {
    size_t i = N[m];
    for (size_t l = 0; l <= m; l++) {
      size_t j = N[l];
      y.fx += Q[i][j] * (y.x[i] * y.x[j] - x.x[i] * x.x[j]);
    }
    for (size_t l = r; l < Q.n; l++) {
      size_t j = N[l];
      y.fx += Q[i][j] * (y.x[i] - x.x[i]) * x.x[j];
    }
  }
}

template <>
void evaluate<evaluation::rfflip_rv>(const ubqp& Q, const solution& x, solution& y, size_t r,
                                     const unique_ptr<size_t[]>& N) {
  y.fx = x.fx;
  for (size_t m = 0; m < r; m++) {
    size_t i = N[m];
    y.fx += (y.x[i] - x.x[i]) * (x.dx[i] + Q[i][i] * (y.x[i] + x.x[i]));
    for (size_t l = 0; l < m; l++) {
      size_t j = N[l];
      y.fx += (y.x[i] - x.x[i]) * (y.x[j] - x.x[j]) * Q[i][j];
    }
  }
}

/**
 * Replaces an incumbent solution with a neighbor solution.
 *
 * @param Q the problem instance.
 * @param eval the evaluation algorithm that was used to evaluate the neighbor.
 * @param x the incumbent solution.
 * @param y the neighbor solution.
 * @param r the number of flip moves between x and y.
 * @param N the permutation vector whose first r indices have been changed between x and y.
 */
template <evaluation eval>
void replace(const ubqp& Q, solution& x, solution& y, size_t r, const unique_ptr<size_t[]>& N) {
  if constexpr (eval == evaluation::f_inc) {
    swap(x.dx, y.dx);
    if (r > 0) {
      size_t i = N[r - 1];
      for (size_t j = 0; j < Q.n; j++)
        if (j != i) x.dx[j] += Q[i][j] * (y.x[i] - x.x[i]);
    }
  } else if constexpr (eval == evaluation::rfflip_rv)
    update_rv(Q, x, y, r, N);
  if constexpr (eval == evaluation::basic)
    swap(x.x, y.x);
  else
    for (size_t k = 0; k < r; k++) x.x[N[k]] = y.x[N[k]];
  x.fx = y.fx;
}

/**
 * Randomizes an incumbent solution.
 *
 * @param Q the problem instance.
 * @param x the incumbent solution to randomize.
 * @param F the maximum value within the components of x.
 * @param bin whether the solution shall be feasible.
 * @param rng the random number generator.
 */
void randomize_incumbent(const ubqp& Q, solution& x, long F, mt19937& rng) {
  x.fx = 0;
  for (size_t i = 0; i < Q.n; i++)
    if ((x.x[i] = (rng() % (F + 1))))
      for (size_t j = 0; j <= i; j++) x.fx += x.x[i] * x.x[j] * Q[i][j];
  assert(check_incumbent_evaluation(Q, x));
  build_rv(Q, x);
  assert(check_rv(Q, x));
}

/**
 * Performs a fractional local search.
 *
 * @param Q the problem instance.
 * @param x the incumbent solution to perform the search on.
 * @param y buffer for neighbor solutions.
 * @param F the maximum value within the components of x and y.
 * @param r the number of flip moves between x and y.
 * @param l_max the maximum amount of iterations without improvements.
 * @param N the permutation vector whose first r indices have been changed between x and y.
 * @param rng the random number generator.
 */
template <evaluation eval>
void ls(const ubqp& Q, solution& x, solution& y, long F, size_t r, size_t l_max,
        unique_ptr<size_t[]>& N, mt19937& rng) {
  for (size_t l = 1; l <= l_max; l++) {
    for (size_t m = 0; m < r; m++) {
      swap(N[m], N[m + (rng() % (Q.n - m))]);
      y.x[N[m]] = rng() % (F + 1);
    }
    evaluate<eval>(Q, x, y, r, N);
    assert(check_evaluation(Q, x, y, r, N));
    if (y < x) {
      replace<eval>(Q, x, y, r, N);
      if constexpr (eval == evaluation::f_inc || eval == evaluation::rfflip_rv)
        assert(check_rv(Q, x));
      assert(check_incumbent_evaluation(Q, x));
      l = 0;
    }
  }
}

/**
 * Represents a relinking method to use in a relinking procedure.
 */
enum struct relinking { delta, random };

/**
 * Performs a relinking procedure.
 *
 * @param rl the relinking method.
 * @param eval the evaluation algorithm to use in the local search.
 * @param Q the problem instance.
 * @param x the initiating solution.
 * @param y the guiding solution.
 * @param w the incumbent output solution.
 * @param z buffer for neighbor solutions.
 * @param F the maximum value within the components of x and y.
 * @param N buffer for permutation vectors.
 * @param rng the random number generator.
 * @return the number of components which are different between x and y.
 */
template <relinking rl, evaluation eval>
size_t relink(const ubqp& Q, const solution& x, const solution& y, solution& w, solution& z, long F,
              unique_ptr<size_t[]>& N, mt19937& rng) {
  (void)F, (void)rng;
  iota(&N[0], &N[Q.n], 0);
  size_t r = 0;
  w.fx = x.fx;
  for (size_t i = 0; i < Q.n; i++) {
    if ((w.x[i] = x.x[i]) != y.x[i]) swap(N[i], N[r++]);
    w.dx[i] = x.dx[i];
  }

  w.fx *= (F * F);
  for (size_t i = 0; i < Q.n; i++) w.x[i] *= F;
  assert(check_incumbent_evaluation(Q, w));
  build_rv(Q, w);
  assert(check_rv(Q, w));

  shuffle(&N[0], &N[r], rng);
  for (size_t k = 0; k < r / 3; k++) z.x[N[k]] = F - w.x[N[k]];
  evaluate<eval>(Q, w, z, r / 3, N);
  assert(check_evaluation(Q, w, z, r / 3, N));

  replace<eval>(Q, w, z, r / 3, N);
  if constexpr (eval == evaluation::f_inc || eval == evaluation::rfflip_rv) assert(check_rv(Q, w));
  assert(check_incumbent_evaluation(Q, w));

  // // if constexpr (rl == relinking::delta)
  // sort(&N[0], &N[r], [&](size_t a, size_t b) {
  //   long fa = a % 2 == 0 ? ((F / 2) * (x.dx[a] + Q[a][a] * (3 * F) / 2))
  //                        : (-(F / 2) * (x.dx[a] + Q[a][a] * F / 2));
  //   long fb = b % 2 == 0 ? ((F / 2) * (x.dx[b] + Q[b][b] * (3 * F) / 2))
  //                        : (-(F / 2) * (x.dx[b] + Q[b][b] * F / 2));
  //   return fa < fb;
  // });
  // // else
  // // for (size_t k = 0; k < r / 3; k++) swap(N[k], N[k + (rng() % (r - k))]);

  swap_ranges(&N[0], &N[r / 3], &N[r - (r / 3)]);
  long f_max = w.fx;
  size_t r_max = 0;
  for (size_t r_curr = 0; r_curr < r / 3; r_curr++) {
    z.x[N[r_curr]] = F - w.x[N[r_curr]];
    // z.x[N[r_curr]] = N[r_curr] % 2 == 0 ? F : 0;
    evaluate<eval>(Q, w, z, r_curr + 1, N);
    assert(check_evaluation(Q, w, z, r_curr + 1, N));
    if (z.fx < f_max) {
      f_max = z.fx;
      r_max = r_curr + 1;
    }
  }

  z.fx = f_max;
  replace<eval>(Q, w, z, r_max, N);
  assert(check_incumbent_evaluation(Q, w));
  if constexpr (eval != evaluation::f_inc && eval != evaluation::rfflip_rv) build_rv(Q, w);
  assert(check_rv(Q, w));

  w.fx /= (F * F);
  for (size_t i = 0; i < Q.n; i++) w.x[i] /= F;
  assert(check_incumbent_evaluation(Q, w));
  build_rv(Q, w);
  assert(check_rv(Q, w));

  return r;
}

using relink_queue = priority_queue<pair<size_t, size_t>, vector<pair<size_t, size_t>>,
                                    greater<pair<size_t, size_t>>>;

bool update_refset(const ubqp& Q, const vector<solution>& ref, size_t size,
                   const unordered_set<size_t>& ref_hashes, const solution& w) {
  (void)ref_hashes;
  return !(any_of(ref.begin(), ref.begin() + size, [&](const solution& x) {
    if (w.fx != x.fx) return false;
    // if (!equal(&w.x[0], &w.x[Q.n], &x.x[0], &x.x[Q.n])) return false;
    return true;
  }));
}

/**
 * Performs a Path Relinking procedure.
 *
 * @param rl the relinking method.
 * @param eval the evaluation algorithm to use in the local search.
 * @param count whether to count the values of `r` throughout relinkings.
 * @param Q the problem instance.
 * @param ref the reference set, where the first is the incumbent.
 * @param ref_hashes set which stores the hashes of the reference set.
 * @param w buffer for relinked solutions.
 * @param z buffer for neighbor solutions.
 * @param F the maximum value within the components of x and y.
 * @param N buffer for permutation vectors.
 * @param L the vector for storing tabu moves.
 * @param K the tabu tenure.
 * @param l_ts_max the maximum amount of iterations for the tabu search procedure.
 * @param l_max the maximum amount of iterations.
 * @param rng the random number generator.
 * @param num_rs the output total number of relinkings.
 * @param sum_rs the sum of the number of values returned by relinkings.
 */
template <relinking rl, evaluation eval, bool count>
void pr(const ubqp& Q, vector<solution>& ref, unordered_set<size_t>& ref_hashes,
        relink_queue& queue, solution& w, solution& z, long F, unique_ptr<size_t[]>& N,
        unique_ptr<size_t[]>& L, size_t K, size_t l_ts_max, size_t l_max, mt19937& rng,
        size_t& num_rs, size_t& sum_rs) {
  for (size_t l = 1; l <= l_max; l++) {
    swap(*min_element(ref.begin(), ref.end()), ref.front());
    ref_hashes.clear();
    ref_hashes.emplace(ref.front().hash(Q));
    for (size_t i = 1; i < ref.size(); i++) {
      size_t hash;
      do {
        // randomize_incumbent(Q, w, F, rng);
        randomize_incumbent(Q, w, 1, rng);
        ts(Q, w, z, F, l_ts_max, L, N, K, rng);
        hash = w.hash(Q);
        // cout << w.fx << endl;
      } while (!update_refset(Q, ref, i, ref_hashes, w));
      ref_hashes.emplace(hash);
      swap(w, ref[i]);
    }
    size_t best = min_element(ref.begin(), ref.end()) - ref.begin();
    size_t worst = max_element(ref.begin(), ref.end()) - ref.begin();
    for (size_t i = 0; i < ref.size(); i++)
      for (size_t j = 0; j < ref.size(); j++)
        if (j != i) queue.emplace(i, j);
    while (queue.size()) {
      size_t i = queue.top().first, j = queue.top().second;
      queue.pop();
      solution& a = ref[i];
      solution& b = ref[j];
      size_t r_max = relink<rl, eval>(Q, a, b, w, z, F, N, rng);
      if constexpr (count) {
        num_rs += 1;
        sum_rs += r_max;
      }
      ts(Q, w, z, F, l_ts_max, L, N, K, rng);
      // cout << ref[best].fx << " " << ref[worst].fx << " ";
      // cout << w.fx << " " << i << ":" << j << " worst " << worst;
      // cout << " " << ref_hashes.size() << " r " << r_max << endl;
      if (w.fx >= ref[worst].fx) continue;
      if (!update_refset(Q, ref, ref.size(), ref_hashes, w)) continue;
      if (w < ref[best]) l = 0;
      ref_hashes.erase(ref[worst].hash(Q));
      ref_hashes.emplace(w.hash(Q));
      swap(w, ref[worst]);
      // cout << "ref";
      // vector<long> v;
      // for (solution& w : ref) v.emplace_back(w.fx);
      // sort(v.begin(), v.end());
      // for (long w : v) cout << " " << w;
      // cout << endl;
      // cout << "improvement" << endl;
      if (worst <= i)
        for (size_t j2 = 0; j2 < (worst == i ? j + 1 : ref.size()); j2++)
          if (j2 != worst) {
            // cout << "emplacing " << worst << ":" << j2 << endl;
            queue.emplace(worst, j2);
          }
      for (size_t i2 = 0; i2 < (worst <= j ? i + 1 : i); i2++)
        if (i2 != worst) {
          // cout << "emplacing " << i2 << ":" << worst << endl;
          queue.emplace(i2, worst);
        }
      best = min_element(ref.begin(), ref.end()) - ref.begin();
      worst = max_element(ref.begin(), ref.end()) - ref.begin();
    }
  }
}

template <evaluation eval>
size_t greedy_start(const ubqp& Q, solution& y, solution& z, long F, unique_ptr<size_t[]>& N,
                    mt19937& rng) {
  randomize_incumbent(Q, y, F, rng);
  assert(check_incumbent_evaluation(Q, y));
  size_t r = 1 + (rng() % (Q.n - 1));
  assert(r > 0 && r <= Q.n);
  for (size_t m = 0; m < r; m++) swap(N[m], N[m + (rng() % (Q.n - m))]);
  for (size_t l = 1; l <= r; l++) {
    for (size_t m = 0; m < r; m++) z.x[N[m]] = rng() % (F + 1);
    evaluate<eval>(Q, y, z, r, N);
    assert(check_evaluation(Q, y, z, r, N));
    if (z < y) {
      replace<eval>(Q, y, z, r, N);
      l = 0;
    }
  }
  for (size_t k = 0; k < r; k++) z.x[N[k]] = y.x[N[k]] < F / 2 ? 0 : F;
  evaluate<eval>(Q, y, z, r, N);
  assert(check_evaluation(Q, y, z, r, N));
  replace<eval>(Q, y, z, r, N);
  if constexpr (eval != evaluation::f_inc && eval != evaluation::rfflip_rv) build_rv(Q, y);
  assert(check_rv(Q, y));
  return r;
}

template <evaluation eval, bool count>
void grasp_ts(const ubqp& Q, solution& x, solution& y, solution& z, long F, unique_ptr<size_t[]>& N,
              unique_ptr<size_t[]>& L, size_t K, size_t l_ts_max, size_t l_max, mt19937& rng,
              size_t& num_rs, size_t& sum_rs) {
  for (size_t l = 1; l <= l_max; l++) {
    size_t r = greedy_start<eval>(Q, y, z, F, N, rng);
    if constexpr (count) {
      num_rs += 1;
      sum_rs += r;
    }
    ts(Q, y, z, F, l_ts_max, L, N, K, rng);
    if (y < x) {
      swap(x, y);
      l = 0;
    }
  }
}

size_t measure(function<void(void)> f) {
  auto t1 = steady_clock::now();
  f();
  auto t2 = steady_clock::now();
  return duration_cast<nanoseconds>(t2 - t1).count();
}

json eval_experiment(const ubqp& Q, json params) {
  mt19937 rng;
  size_t r = params["r"];
  long F = params["F"];
  size_t l_max = params["l_max"];
  evaluation eval = static_cast<evaluation>(params["eval"]);
  unique_ptr<size_t[]> N(new size_t[Q.n]);
  iota(&N[0], &N[Q.n], 0);
  solution x(Q), y(Q);
  long average = 0, maximum = numeric_limits<long>::min(), minimum = numeric_limits<long>::max();
  for (size_t i = 0; i < l_max; i++) {
    randomize_incumbent(Q, x, F, rng);
    for (size_t m = 0; m < r; m++) {
      swap(N[m], N[m + (rng() % (Q.n - m))]);
      y.x[N[m]] = rng() % (F + 1);
    }
    long dt = 0;
    switch (eval) {
      case evaluation::basic:
        dt = measure([&]() { evaluate<evaluation::basic>(Q, x, y, r, N); });
        break;
      case evaluation::f_inc:
        dt = measure([&]() { evaluate<evaluation::f_inc>(Q, x, y, r, N); });
        break;
      case evaluation::rfflip:
        dt = measure([&]() { evaluate<evaluation::rfflip>(Q, x, y, r, N); });
        break;
      case evaluation::rfflip_rv:
        dt = measure([&]() { evaluate<evaluation::rfflip_rv>(Q, x, y, r, N); });
        break;
    }
    average += dt;
    if (dt > maximum) maximum = dt;
    if (dt < minimum) minimum = dt;
  }
  average /= l_max;
  return {{"avg", average}, {"max", maximum}, {"min", minimum}};
}

json ls_experiment(const ubqp& Q, json params) {
  mt19937 rng;
  size_t r = params["r"];
  long F = params["F"];
  size_t l_max = params["l_max"];
  evaluation eval = static_cast<evaluation>(params["eval"]);
  unique_ptr<size_t[]> N(new size_t[Q.n]);
  iota(&N[0], &N[Q.n], 0);
  solution x(Q), y(Q);
  randomize_incumbent(Q, x, F, rng);
  long dt = 0;
  switch (eval) {
    case evaluation::basic:
      dt = measure([&]() { ls<evaluation::basic>(Q, x, y, F, r, l_max, N, rng); });
      break;
    case evaluation::f_inc:
      dt = measure([&]() { ls<evaluation::f_inc>(Q, x, y, F, r, l_max, N, rng); });
      break;
    case evaluation::rfflip:
      dt = measure([&]() { ls<evaluation::rfflip>(Q, x, y, F, r, l_max, N, rng); });
      break;
    case evaluation::rfflip_rv:
      dt = measure([&]() { ls<evaluation::rfflip_rv>(Q, x, y, F, r, l_max, N, rng); });
      break;
  }
  return {{"dt", dt}, {"fx", x.fx / (F * F)}};
}

json pr_experiment(const ubqp& Q, json params) {
  mt19937 rng;
  long F = params["F"];
  size_t K = params["K"];
  size_t ref_size = params["ref_size"];
  size_t l_ts_max = params["l_ts_max"];
  size_t l_max = params["l_max"];
  bool count = params["count"];
  relinking rm = static_cast<relinking>(params["relink"]);
  evaluation eval = static_cast<evaluation>(params["eval"]);
  size_t num_rs = 0, sum_rs = 0;
  unique_ptr<size_t[]> N(new size_t[Q.n]), L(new size_t[Q.n]);
  iota(&N[0], &N[Q.n], 0);
  solution w(Q), z(Q);
  vector<solution> ref;
  unordered_set<size_t> ref_hashes;
  relink_queue queue;
  ref_hashes.reserve(ref_size);
  for (size_t i = 0; i < ref_size; i++) ref.emplace_back(Q);
  randomize_incumbent(Q, ref.front(), 1, rng);
  // randomize_incumbent(Q, ref.front(), F, rng);
  ts(Q, ref.front(), z, F, l_ts_max, L, N, K, rng);
  if (!count) {
    long dt = 0;
    switch (rm) {
      case relinking::delta:
        switch (eval) {
          case evaluation::basic:
            dt = measure([&]() {
              pr<relinking::delta, evaluation::basic, false>(Q, ref, ref_hashes, queue, w, z, F, N,
                                                             L, K, l_ts_max, l_max, rng, num_rs,
                                                             sum_rs);
            });
            break;
          case evaluation::f_inc:
            dt = measure([&]() {
              pr<relinking::delta, evaluation::f_inc, false>(Q, ref, ref_hashes, queue, w, z, F, N,
                                                             L, K, l_ts_max, l_max, rng, num_rs,
                                                             sum_rs);
            });
            break;
          case evaluation::rfflip:
            dt = measure([&]() {
              pr<relinking::delta, evaluation::rfflip, false>(Q, ref, ref_hashes, queue, w, z, F, N,
                                                              L, K, l_ts_max, l_max, rng, num_rs,
                                                              sum_rs);
            });
            break;
          case evaluation::rfflip_rv:
            dt = measure([&]() {
              pr<relinking::delta, evaluation::rfflip_rv, false>(Q, ref, ref_hashes, queue, w, z, F,
                                                                 N, L, K, l_ts_max, l_max, rng,
                                                                 num_rs, sum_rs);
            });
            break;
        }
        break;
      case relinking::random:
        switch (eval) {
          case evaluation::basic:
            dt = measure([&]() {
              pr<relinking::random, evaluation::basic, false>(Q, ref, ref_hashes, queue, w, z, F, N,
                                                              L, K, l_ts_max, l_max, rng, num_rs,
                                                              sum_rs);
            });
            break;
          case evaluation::f_inc:
            dt = measure([&]() {
              pr<relinking::random, evaluation::f_inc, false>(Q, ref, ref_hashes, queue, w, z, F, N,
                                                              L, K, l_ts_max, l_max, rng, num_rs,
                                                              sum_rs);
            });
            break;
          case evaluation::rfflip:
            dt = measure([&]() {
              pr<relinking::random, evaluation::rfflip, false>(Q, ref, ref_hashes, queue, w, z, F,
                                                               N, L, K, l_ts_max, l_max, rng,
                                                               num_rs, sum_rs);
            });
            break;
          case evaluation::rfflip_rv:
            dt = measure([&]() {
              pr<relinking::random, evaluation::rfflip_rv, false>(Q, ref, ref_hashes, queue, w, z,
                                                                  F, N, L, K, l_ts_max, l_max, rng,
                                                                  num_rs, sum_rs);
            });
            break;
        }
        break;
    }
    return {{"fx", min_element(ref.begin(), ref.end())->fx / (F * F)}, {"dt", dt}};
  } else {
    switch (rm) {
      case relinking::delta:
        switch (eval) {
          case evaluation::basic:
            pr<relinking::delta, evaluation::basic, true>(Q, ref, ref_hashes, queue, w, z, F, N, L,
                                                          K, l_ts_max, l_max, rng, num_rs, sum_rs);
            break;
          case evaluation::f_inc:
            pr<relinking::delta, evaluation::f_inc, true>(Q, ref, ref_hashes, queue, w, z, F, N, L,
                                                          K, l_ts_max, l_max, rng, num_rs, sum_rs);
            break;
          case evaluation::rfflip:
            pr<relinking::delta, evaluation::rfflip, true>(Q, ref, ref_hashes, queue, w, z, F, N, L,
                                                           K, l_ts_max, l_max, rng, num_rs, sum_rs);
            break;
          case evaluation::rfflip_rv:
            pr<relinking::delta, evaluation::rfflip_rv, true>(
                Q, ref, ref_hashes, queue, w, z, F, N, L, K, l_ts_max, l_max, rng, num_rs, sum_rs);
            break;
        }
        break;
      case relinking::random:
        switch (eval) {
          case evaluation::basic:
            pr<relinking::random, evaluation::basic, true>(Q, ref, ref_hashes, queue, w, z, F, N, L,
                                                           K, l_ts_max, l_max, rng, num_rs, sum_rs);
            break;
          case evaluation::f_inc:
            pr<relinking::random, evaluation::f_inc, true>(Q, ref, ref_hashes, queue, w, z, F, N, L,
                                                           K, l_ts_max, l_max, rng, num_rs, sum_rs);
            break;
          case evaluation::rfflip:
            pr<relinking::random, evaluation::rfflip, true>(
                Q, ref, ref_hashes, queue, w, z, F, N, L, K, l_ts_max, l_max, rng, num_rs, sum_rs);
            break;
          case evaluation::rfflip_rv:
            pr<relinking::random, evaluation::rfflip_rv, true>(
                Q, ref, ref_hashes, queue, w, z, F, N, L, K, l_ts_max, l_max, rng, num_rs, sum_rs);
            break;
        }
        break;
    }
    return {{"num_rs", num_rs}, {"sum_rs", sum_rs}};
  }
}

int main(int argc, const char* argv[]) {
  (void)argc, (void)argv;

  // long fx_avg = 0;
  // for (string instance : {
  //          "G1",        "G2",        "G3",        "G4",        "G5",         "G6",
  //          "G7",        "G8",        "G9",        "G10",       "G11",        "G12",
  //          "G13",       "G14",       "G15",       "G16",       "G17",        "G18",
  //          "G19",       "G20",       "G21",       "G22",       "G23",        "G24",
  //          "G25",       "G25",       "G26",       "G27",       "G28",        "G29",
  //          "G30",       "G31",       "G32",       "G33",       "G34",        "G35",
  //          "G36",       "G37",       "G38",       "G39",       "G40",        "G41",
  //          "G42",       "bqp2500.1", "bqp2500.2", "bqp2500.3", "bqp2500.4",  "bqp2500.5",
  //          "bqp2500.6", "bqp2500.7", "bqp2500.8", "bqp2500.9", "bqp2500.10",
  //      }) {
  //   auto Q = ubqp::load(instance);
  //   mt19937 rng;
  //   solution inc(Q), cur(Q);
  //   auto N = make_unique<size_t[]>(Q.n);
  //   for (size_t i = 1; i <= 2; i++) {
  //     randomize_incumbent(Q, inc, 1, rng);
  //     ts(Q, inc, cur, 0, 0, N, N, 0, rng);
  //     cout << instance << " " << inc.fx << endl;
  //     fx_avg += inc.fx;
  //   }
  // }
  // cout << "average " << fx_avg / 106 << endl;
  // return 0;

  if (argc < 3) return 1;
  const string results_filename(argv[1]);
  json results;
  mutex results_mutex;
  {
    ifstream file(results_filename, ifstream::in);
    if (file) file >> results;
  }

  auto run = [&](function<json(const ubqp&, json)> exp, const ubqp& Q, json params) {
    {
      lock_guard<mutex> lock(results_mutex);
      for (auto& j : results)
        if (j["params"] == params) return;
      cout << "<<< " << params << endl;
    }
    json result = exp(Q, params);
    {
      lock_guard<mutex> lock(results_mutex);
      cout << ">>> " << result << endl;
      results.push_back({{"params", params}, {"result", result}});
      ofstream file(results_filename, ifstream::out);
      file << setw(1) << results;
    }
  };

  const string experiment(argv[2]);

  if (experiment == "eval")
    for (string instance : {"bqp1000.1", "bqp2500.1", "G43", "G22"}) {
      ubqp Q = ubqp::load(instance);
      size_t r_step = max(Q.n / 100, size_t(1));
      for (evaluation eval : {
               evaluation::basic,
               evaluation::f_inc,
               evaluation::rfflip,
               evaluation::rfflip_rv,
           })
        for (long F : {10, 100})
          for (size_t r = r_step; r <= Q.n; r += r_step)
            run(eval_experiment, Q,
                {{"exp", "eval"},
                 {"instance", instance},
                 {"n", Q.n},
                 {"eval", static_cast<int>(eval)},
                 {"r", r},
                 {"F", F},
                 {"l_max", 100}});
    }

  if (experiment == "ls")
    for (string instance : {"bqp1000.1", "bqp2500.1", "G43", "G22"}) {
      ubqp Q = ubqp::load(instance);
      size_t r_step = max(Q.n / 100, size_t(1));
      for (evaluation eval : {
               evaluation::basic,
               evaluation::f_inc,
               evaluation::rfflip,
               evaluation::rfflip_rv,
           })
        for (long F : {10, 100})
          for (size_t r = r_step; r <= Q.n; r += r_step)
            run(ls_experiment, Q,
                {{"exp", "ls"},
                 {"instance", instance},
                 {"n", Q.n},
                 {"eval", static_cast<int>(eval)},
                 {"r", r},
                 {"F", F},
                 {"l_max", Q.n}});
    }

  vector<future<void>> futures;

  if (experiment == "pr" || experiment == "pr_count") {
    bool count = experiment == "pr_count";
    for (string instance : {
             "bqp250.1",  "bqp250.2",   "bqp250.3",   "bqp250.4",  "bqp250.5",  "bqp250.6",
             "bqp250.7",  "bqp250.8",   "bqp250.9",   "bqp250.10", "bqp500.1",  "bqp500.2",
             "bqp500.3",  "bqp500.4",   "bqp500.5",   "bqp500.6",  "bqp500.7",  "bqp500.8",
             "bqp500.9",  "bqp500.10",  "G1",         "G2",        "G3",        "G4",
             "G5",        "G6",         "G7",         "G8",        "G9",        "G10",
             "G11",       "G12",        "G13",        "G14",       "G15",       "G16",
             "G17",       "G18",        "G19",        "G20",       "G21",       "bqp1000.1",
             "bqp1000.2", "bqp1000.3",  "bqp1000.4",  "bqp1000.5", "bqp1000.6", "bqp1000.7",
             "bqp1000.8", "bqp1000.9",  "bqp1000.10", "G43",       "G44",       "G45",
             "G46",       "G47",        "G51",        "G52",       "G53",       "G54",
             "G22",       "G23",        "G24",        "G25",       "G25",       "G26",
             "G27",       "G28",        "G29",        "G30",       "G31",       "G32",
             "G33",       "G34",        "G35",        "G36",       "G37",       "G38",
             "G39",       "G40",        "G41",        "G42",       "bqp2500.1", "bqp2500.2",
             "bqp2500.3", "bqp2500.4",  "bqp2500.5",  "bqp2500.6", "bqp2500.7", "bqp2500.8",
             "bqp2500.9", "bqp2500.10",
         }) {
      auto Q = make_shared<ubqp>(move(ubqp::load(instance)));
      for (relinking rl : {relinking::random}) {
        for (evaluation eval : {
                 evaluation::basic,
                 evaluation::f_inc,
                 evaluation::rfflip,
                 evaluation::rfflip_rv,
             })
          if (!count || eval == evaluation::rfflip_rv)
            futures.emplace_back(async([=]() {
              run(pr_experiment, *Q,
                  {{"exp", "pr"},
                   {"instance", instance},
                   {"n", Q->n},
                   {"count", count},
                   {"eval", static_cast<int>(eval)},
                   {"relink", static_cast<int>(rl)},
                   {"F", 10},
                   {"K", Q->n / 10},
                   {"ref_size", 10},
                   {"l_ts_max", 10000},
                   {"l_max", 10}});
            }));
        if (!count) {
          for (auto& f : futures) f.wait();
          futures.clear();
        }
      }
    }
    for (auto& f : futures) f.wait();
    futures.clear();
  }

  return 0;
}
