/* inc/LHSPoint.hpp
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


/* Generates a random point by the latin hyperbube sampling technique 
* 
* Given hyper-rectangle of dimension DIM and N requested number of 
* sample points, generate a DIM x N matrix with each row being a 
** random permutation of the row vector [1, 2, ..., N] 
* 
* 
* The random integers are generated by the algorithm 
* implement by John Burkardt
*
*  Scholarly reference:
*
*    Brian Beachkofski, Ramana Grandhi,
*    Improved Distributed Hypercube Sampling,
*    American Institute of Aeronautics and Astronautics Paper 2002-1274.
*
*/

#ifndef LHSPOINT_HPP
#define LHSPOINT_HPP

#include <random>
#include <vector>

#include "ihs.hpp"
#include "DynProg.hpp"

namespace DTDP
{
class LHSPoint
{
public:
  // using D = double;
  // using Index = int;
  using unif_dist = std::uniform_int_distribution<Index>;

  LHSPoint() = default;
  LHSPoint(const LHSPoint &) = default;
  LHSPoint &operator=(const LHSPoint &) = default;
  LHSPoint(LHSPoint &&) = default;
  LHSPoint &operator=(LHSPoint &&) = default;
  ~LHSPoint() = default;

  /* pt_dim: dimension of the point, = dim of hypercube rectangle
  pt_count: number of points to be generated
  duplication_fac: duplication factor, d >= 1
  */
  LHSPoint(Index pt_dim,
           Index pt_count,
           Index duplication_fac = 5,
           Index seed = 1543);

  /* if a point is provided around which to sample a nearby point, 
  * one can provide the following pertubation factors
  * 
  */

  constexpr static D default_del_u{0.0075}; /* bounds are set to +/- it if pt(i) = 0 */
  constexpr static D default_del_z{0.05};   /* and, otherwise, to pt(i) +/- del_z*pt(i) */

  // constexpr static D default_del_u{0.075}; /* bounds are set to +/- it if pt(i) = 0 */
  // constexpr static D default_del_z{0.25};   /* and, otherwise, to pt(i) +/- del_z*pt(i) */

  /* returns near_pt, a perturbation point near pt 
  * given perturbation factors del_u and del_z
  * as determinants of lower and upper bounds of the variable
  */
  void get_nearby_pt(const D *pt,
                     D *near_pt,
                     Index pt_dim,
                     Index j_lhs_pt, /* point ID in the lhs */
                     D del_u = default_del_u,
                     D del_z = default_del_z);

  /* returns near_pt, a perturbation point near pt 
  * given lower and upper bounds of the variable
  */
  void get_nearby_pt(D *near_pt,
                     const D *lb, /* lower bound */
                     const D *ub, /* upper bound */
                     Index pt_dim,
                     Index j_lhs_pt /* point ID in the lhs */);

  Index get_pt_count() const { return pt_count; }

protected:
  void get_nearby_pt(D *near_pt,
                     Index pt_dim,
                     D lb,
                     D ub,
                     Index i_pt_elem, /* ID of element in pt */
                     Index j_lhs_pt /* point ID in the lhs */)
  {
    near_pt[i_pt_elem] = lb + rand(i_pt_elem, j_lhs_pt, pt_dim) * (ub - lb);
  }

  D rand(Index i_pt_elem, Index j_lhs_pt, Index pt_dim);

private:
  std::vector<Index> lhs_x;
  Index pt_count;
  std::mt19937 mtrng;
  unif_dist udist;
};
} // namespace DTDP

#endif
