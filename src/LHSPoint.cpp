
/* src/LHSPoint.cpp
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


#include <algorithm>

#include "LHSPoint.hpp"
#include "Exception.hpp"

namespace DTDP
{
LHSPoint::LHSPoint(Index pt_dim, Index pt_count, Index dup_fac, Index seed)
    : lhs_x{std::vector<Index>(pt_dim * pt_count, 0)},
      pt_count{pt_count},
      mtrng{std::mt19937(seed)},
      udist{unif_dist(1, pt_count)}
{

  if (pt_dim < 1)
    Exception::error(__FILE__, __LINE__, std::string{"dim < 1."});

  if (pt_count < 1)
    Exception::error(__FILE__, __LINE__, std::string{"n < 1."});

  if (seed < 0)
    Exception::error(__FILE__, __LINE__, std::string{"seed < 0."});

  Index xsize = lhs_x.size();
  Index *x = new Index[xsize];

  x = ihs(pt_dim, pt_count, dup_fac, seed);

  std::copy(&x[0], &x[0] + xsize, std::begin(lhs_x));

  delete[] x;
}

void LHSPoint::get_nearby_pt(const D *pt,
                             D *near_pt,
                             Index pt_dim,
                             Index j_lhs_pt, /* point ID in the lhs */
                             D del_u,
                             D del_z)
{
  // returns near_pt, a perturbation point near pt
  D lb, ub;
  for (Index i_pt_elem{0}; i_pt_elem < pt_dim; ++i_pt_elem)
  {
    // get bounds
    if (pt[i_pt_elem] < 0.0)
    {
      lb = (1.0 + del_z) * pt[i_pt_elem];
      ub = (1.0 - del_z) * pt[i_pt_elem];
    }
    else if (pt[i_pt_elem] > 0.0)
    {
      lb = (1.0 - del_z) * pt[i_pt_elem];
      ub = (1.0 + del_z) * pt[i_pt_elem];
    }
    else
    {
      lb = (0.0 - del_u);
      ub = del_u;
    }

    // get_bounds(lb, ub, pt, i_pt_elem, del_u, del_z);
    get_nearby_pt(&near_pt[0], pt_dim, lb, ub, i_pt_elem, j_lhs_pt);
  }
}

void LHSPoint::get_nearby_pt(D *near_pt,
                             const D *lb, /* lower bound */
                             const D *ub, /* upper bound */
                             Index pt_dim,
                             Index j_lhs_pt /* point ID in the lhs */)
{
  // returns near_pt, a perturbation point near pt
  for (Index i{0}; i < pt_dim; ++i)
  {
    get_nearby_pt(&near_pt[0],
                  pt_dim,
                  lb[i],
                  ub[i],
                  i,
                  j_lhs_pt);
  }
}

D LHSPoint::rand(Index i_pt_elem, Index j_lhs_pt, Index pt_dim)
{
  const Index rx{lhs_x[i_pt_elem + j_lhs_pt * pt_dim]};
  const Index rx1{udist(mtrng)};
  return static_cast<D>(rx - rx1) / static_cast<D>(pt_count);
}
} // namespace DTDP