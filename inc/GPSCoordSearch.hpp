/* inc/GPSCoordSearch.hpp
*  - Generalized pattern search (GPS)
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

#ifndef GPSCoordSearch_HPP
#define GPSCoordSearch_HPP

#include <random>

#include "CoordSearch.hpp"
#include "LHSPoint.hpp"

namespace DTDP
{
class GPSCoordSearch : virtual public CoordSearch,
                       virtual public LHSPoint
{
public:
  using D = CoordSearch::D;
  using Index = CoordSearch::Index;

  GPSCoordSearch() = default;
  GPSCoordSearch(const GPSCoordSearch &) = default;
  GPSCoordSearch &operator=(const GPSCoordSearch &) = default;
  virtual ~GPSCoordSearch() = default;

  explicit GPSCoordSearch(const D *x0,
                          Index n,
                          Index pattern_type,
                          D line_search_coeff,
                          Index lhs_pt_count, /* number of points to generate if latin hypercube sampling */
                          Index search_method);

  explicit GPSCoordSearch(const D *x0,
                          Index n,
                          D delta_init,
                          D delta_tol,
                          Index it_max,
                          Index fc_max,
                          D ecoeff,
                          D ccoeff,
                          Index pattern_type,
                          D line_search_coeff,
                          Index lhs_pt_count,
                          Index search_method);

  /* virtual methods */
  virtual void init_fval() override;
  virtual void search() override;
  virtual void lhs_search();
  // virtual void update_pattern() override { return; /* DON'T pattern.sort();*/ };
  // virtual void set_bounds();

  void line_search();
  bool lhs_search(Index idx_row, Index &fcalls);
  bool accept_search_pt();

  D get_ls_coeff() const { return ls_coeff; }
  Index get_search_method() const { return search_method; }

  enum class SearchMethod
  {
    LINE_SEARCH = 0,
    LHS_SEARCH, /* latin hypercube sampling */
    NONE
  };

  void verify_search_method();
  void print_lsh_x(std::ostream &os = std::cout) const;

private:
  D ls_coeff; /* if unequal zero, a line search method will be used to select the next iterate */
  Index search_method;
};
} // namespace DTDP

#endif