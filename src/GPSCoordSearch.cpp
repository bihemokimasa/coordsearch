/* src/GPSCoordSearch.cpp
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


#include <iomanip>

#include "GPSCoordSearch.hpp"
#include "Exception.hpp"
#include "ihs.hpp"

namespace DTDP
{
GPSCoordSearch::
    GPSCoordSearch(const D *x0,
                   Index n,
                   Index pattern_type,
                   D line_search_coeff,
                   Index lhs_pt_count,
                   Index search_meth)
    : CoordSearch(&x0[0], n),
      LHSPoint{(search_meth ==
                static_cast<Index>(SearchMethod::LHS_SEARCH))
                   ? LHSPoint(n, lhs_pt_count)
                   : LHSPoint()},
      ls_coeff{line_search_coeff},
      search_method{search_meth}
{
  verify_search_method();
  const D default_ecoeff{2.0};
  set_ecoeff(default_ecoeff); 
  set_pattern(pattern_type);
}

GPSCoordSearch::
    GPSCoordSearch(const D *x0,
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
                   Index search_meth)
    : CoordSearch(&x0[0],
                  n,
                  delta_init,
                  delta_tol,
                  it_max,
                  fc_max,
                  ecoeff,
                  ccoeff),
      LHSPoint{(search_meth ==
                static_cast<Index>(SearchMethod::LHS_SEARCH))
                   ? LHSPoint(n, lhs_pt_count)
                   : LHSPoint()},
      ls_coeff{line_search_coeff},
      search_method{search_meth}
{
  verify_search_method();
  set_pattern(pattern_type);
}

void GPSCoordSearch::verify_search_method()
{
  if ((search_method != static_cast<Index>(SearchMethod::LINE_SEARCH) &&
       search_method != static_cast<Index>(SearchMethod::LHS_SEARCH)) &&
      search_method != static_cast<Index>(SearchMethod::NONE))
    Exception::error(__FILE__, __LINE__,
                     std::string{"choose either 0 (line search)"} +
                         std::string{" or 1 (latin hs search) for search method."});
}

void GPSCoordSearch::init_fval()
{
  /* initialize fx */
  if (search_method == static_cast<Index>(SearchMethod::LHS_SEARCH))
  {
    CoordSearch::init_fval(); // evaluate the first point
    lhs_search();
  }
  else
  {
    CoordSearch::init_fval();
  }
}

void GPSCoordSearch::search()
{
  print_progress();
  init_fval();
  print_progress();
  while (!stop())
  {
    /* try search step */
    if (get_iters() > 0)
    {
      if (search_method == static_cast<Index>(SearchMethod::LINE_SEARCH))
      {
        line_search();
      }
      else if (search_method == static_cast<Index>(SearchMethod::LHS_SEARCH))
      {
        lhs_search();
      }
    }

    if (stop())
    {
      break;
    }

    /* search step wasn't succesfful, try poll step */
    poll();

    update_pattern();
    inc_iters();
    print_progress();
  }
  print_result();
}

void GPSCoordSearch::line_search()
{
  for (Index i{0}; i < get_n(); ++i)
  {
    const D val{get_pt()[i] - get_ptbefore()[i]};
    set_ptnew(i, get_pt()[i] + ls_coeff * val);
  }

  set_fxnew(fval(get_ptnew()));
  inc_fcalls();
  accept_search_pt();
}

bool GPSCoordSearch::
    lhs_search(Index idx_row,
               Index &fcalls)
{
  /*
  * return true to signal termination
  */

  fcalls = 0;

  if (get_fcalls() >= get_fcalls_max() ||
      idx_row >= get_pt_count()) /* to allow for parallelization */
  {
    return true;
  }

  LHSPoint::get_nearby_pt(get_ptnew().data(),
                          get_ptnew().data(),
                          get_n(), idx_row);

  set_fxnew(fval(get_ptnew()));

  fcalls++;

  return false;
}

void GPSCoordSearch::lhs_search()
{
  /* search based on the latin hypercube sampling */

  Index i_row{0};
  Index fcalls{0};
  while (1)
  {
    if (lhs_search(i_row++,
                   fcalls))
    {
      break;
    }

    inc_fcalls(fcalls);
    accept_search_pt();
  }
}

bool GPSCoordSearch::accept_search_pt()
{
  if (get_fxnew() < get_fx())
  {
    set_fx(get_fxnew());
    set_pt(get_ptnew());
  }
}

} // namespace DTDP