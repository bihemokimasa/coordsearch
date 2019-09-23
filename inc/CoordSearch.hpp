/* inc/CoordSearch.hpp
* 
* 
* Copyright (C) 2019 Bihemo Kimasa
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.

* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

/*
* References:
* Virginia Torczon. 1997. 
  On the Convergence of Pattern Search Algorithms. 
  SIAM J. on Optimization 7, 1 (January 1997), 1-25.
*
*
* Audet, Charles & Hare, Warren. (2017). 
* Derivative-Free and Blackbox Optimization. 
*/

#ifndef COORDSEARCH_HPP
#define COORDSEARCH_HPP

#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>
#include <fstream>

#include "testfunc.hpp"

namespace DTDP
{
class CoordSearch
{
public:
  using D = double;
  using Index = int;

  CoordSearch() = default;
  CoordSearch(const CoordSearch &) = default;
  CoordSearch &operator=(const CoordSearch &) = default;
  virtual ~CoordSearch() = default;

  explicit CoordSearch(const D *x0, Index n);

  explicit CoordSearch(const D *x0, Index n,
                       D delta_init, D delta_tol,
                       Index it_max, Index fc_max,
                       D ecoeff, D ccoeff);

  /* virtual methods */
  virtual void init_fval()
  { /* initialize fx */
    set_fx(fval(get_pt()));
    inc_fcalls();
  }

  virtual void search();
  virtual void new_xiter(const std::vector<D> &xbefore, Index i_col, std::vector<D> &xnew);
  virtual void update_pattern() { pattern.sort(); };
  virtual D fval(const D *x, Index n) const { return Rosebrock(&x[0], n); }
  virtual Index pattern_size() const { return pattern.get_col_count(); }
  virtual void inc_fcalls(Index step = 1) { fc += step; }
  virtual void poll();
  virtual bool stop();

  virtual void print_pattern(std::ostream &os = std::cout) const
  {
    os << pattern << std::endl;
  }

  virtual void print_progress(std::ostream &os = std::cout);
  virtual void print_result(std::ostream &os = std::cout);

  void update_slen();

  Index get_row_count(Index i_col) const { return pattern.get_row_count(i_col); }
  void set_pattern_fx(Index i_col, D val) { pattern.pat[i_col].first = val; }

  D fval(const std::vector<D> &x) const { return fval(x.data(), x.size()); }

  bool poll(Index i, Index &fcalls);
  bool accept_ptnew();

  /* attribute access */
  const std::vector<D> &get_pt() const { return pt; }
  std::vector<D> &get_pt() { return pt; }

  D get_fx() const { return fx; }

  D get_slen() const { return sl; }
  Index get_iters() const { return it; }
  Index get_fcalls() const { return fc; }
  Index exit_flag() const { return xt_f; }
  Index get_n() const { return pt.size(); }

  /* return parameters */
  D get_slen_init() const { return par.sl0; }
  D get_slen_min() const { return par.smin; }
  Index get_it_max() const { return par.it_max; }
  Index get_fcalls_max() const { return par.fc_max; }

  /* setters */
  void set_pattern(Index ptype) { pattern = Pattern(get_n(), ptype); }
  void set_pt(const D *x, Index n) { std::copy(&x[0], &x[0] + n, std::begin(pt)); }
  void set_pt(const std::vector<D> &x) { pt = x; }
  void set_fx(D val) { fx = val; }
  void set_slen(D val) { sl = val; };
  void set_iters(Index val) { it = val; }
  void set_fcalls(Index val) { fc = val; }
  void set_exit_flag(Index val) { xt_f = val; }

  void inc_iters(Index step = 1) { it += step; }

  void set_init_slen(D val) { par.sl0 = val; }
  void set_min_slen(D val) { par.smin = val; }
  void set_ecoeff(D val) { par.eco = val; }
  void set_ccoeff(D val) { par.cco = val; }
  void set_max_iters(Index val) { par.it_max = val; }
  void set_max_fcalls(Index val) { par.fc_max = val; }

  bool dont_stop();

  std::vector<D> copy(const D *x, Index n) const;
  void copy(const D *x, Index n, std::vector<D> &dest);
  void copy(const std::vector<D> &src, std::vector<D> &dest);

  /* static variables or functions */
  static D get_infinity() { return std::numeric_limits<double>::infinity(); }
  static Index check_range_and_get_idx(Index i,
                                       Index low = 0,
                                       Index high = std::numeric_limits<int>::max());

  struct DefaultParameters
  {

    DefaultParameters() = default;
    DefaultParameters(const DefaultParameters &) = default;
    DefaultParameters &operator=(const DefaultParameters &) = default;
    ~DefaultParameters() = default;

    constexpr static D sl0{1.0};
    constexpr static D smin{1e-6};
    constexpr static D eco{1.0};
    constexpr static D cc0{0.5};
  } def;

  struct Pattern
  {
    /* a pattern is the matrix P = [I -I] where 
    I is the identiy matrix of dimension n x n
    such pattern is here argumented by the 
    point function value so that its columns 
    may be ranked

    This pattern contains the set of directions when polling

    Don't store the zeros by having 
    each unit vector identified by its row and a value in {1, -1} 
    */

    Pattern() = default;
    Pattern(const Pattern &) = default;
    Pattern &operator=(const Pattern &) = default;
    ~Pattern() = default;

    template <class T>
    using vec = std::vector<T>;

    typedef std::pair<vec<Index>, vec<D>> unity_coord; /* unit coordinate vector, gives search direction */
    vec<std::pair<D, unity_coord>> pat;
    Index poll_type;

    enum class Type
    {
      /* polling type */
      GPS2N = 0, /* P = [I -I]*/
      GPSNp1,    /* P = [I e] where e = (-1, -1, ...,) is an n-column vector of negative ones */
    };

    Pattern(Index n, Index ptype = static_cast<Index>(Type::GPS2N));
    static Index check_ptype_and_return(Index ptype);

    struct less_than
    {
      bool operator()(const std::pair<D, unity_coord> &a,
                      const std::pair<D, unity_coord> &b);
    };

    void sort() { std::sort(std::begin(pat), std::end(pat), less_than()); }
    friend std::ostream &operator<<(std::ostream &os, const Pattern &p);

    Index get_col_count() const { return pat.size(); }
    Index get_row_count(Index i_col) const { return pat[i_col].second.first.size(); }
  };

protected:
  const std::vector<D> &get_ptbefore() const { return ptbefore; }
  std::vector<D> &get_ptbefore() { return ptbefore; }

  D get_fxbefore() const { return fxbefore; }

  void set_ptbefore(const std::vector<D> &x) { ptbefore = x; }
  void set_fxbefore(D val) { fxbefore = val; }

  const std::vector<D> &get_ptnew() const { return ptnew; }
  std::vector<D> &get_ptnew() { return ptnew; }

  D get_fxnew() const { return fxnew; }

  void set_ptnew(const std::vector<D> &x) { ptnew = x; }
  void set_ptnew(const D *x, Index n) { std::copy(&x[0], &x[0] + n, std::begin(ptnew)); }
  void set_ptnew(Index i, D val) { ptnew[i] = val; }
  void set_fxnew(D val) { fxnew = val; }

  void set_fxptbefore();
  void update_fxptbefore();
  void update_fxptnew();
  void update();

private:
  Pattern pattern;
  std::vector<D> pt;              /* function input point vector */
  std::vector<D> ptnew, ptbefore; /* workspace variables: function input point vector */
  D fx;                           /* point function value */
  D fxnew, fxbefore;              /* workspace variables: point function value */
  D sl;                           /* current step length */
  Index it;                       /* current iteration count */
  Index fc;                       /* function calls counter */
  Index xt_f;                     /* exit flag */

  struct Parameters
  {
    D sl0;        /* initial step length */
    D smin;       /* minimum value for the step length */
    D eco;        /* steplegth adjustment parameter: expansion coefficient */
    D cco;        /* steplegth adjustment parameter: contraction coeffient */
    Index it_max; /* maximum iteration count */
    Index fc_max; /* maximum function call count */

    Parameters() = default;
    Parameters(const Parameters &) = default;
    Parameters &operator=(const Parameters &) = default;
    ~Parameters() = default;

    Parameters(D s_init, D s_eps,
               Index it, Index fc,
               D ecoeff, D ccoeff)
        : sl0{s_init}, smin{s_eps},
          eco{ecoeff}, cco{ccoeff},
          it_max{it}, fc_max{fc} {}
  } par;
};

} // namespace DTDP

#endif