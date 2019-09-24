/* src/CoordSearch.cpp
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

#include <random>
#include <iomanip>

#include "CoordSearch.hpp"
#include "Exception.hpp"

namespace DTDP
{

CoordSearch::CoordSearch(const D *x0, Index n)
    : CoordSearch(&x0[0], check_range_and_get_idx(n, 1),
                  def.sl0,
                  def.smin,
                  n * 20000,
                  n * 100,
                  def.eco,
                  def.cc0) {}

CoordSearch::CoordSearch(const D *x0,
                         Index n,
                         D s_init,
                         D s_tol,
                         Index it_max,
                         Index fc_max,
                         D ecoeff,
                         D ccoeff)
    : pattern(n),
      pt{copy(&x0[0], check_range_and_get_idx(n, 1))},
      ptnew{std::vector<D>(n, get_infinity())},
      ptbefore{std::vector<D>(n, get_infinity())},
      fx{get_infinity()},
      fxnew{get_infinity()},
      fxbefore{get_infinity()},
      sl{s_init},
      it{0},
      fc{0},
      xt_f{-1},
      par{Parameters(s_init,
                     s_tol,
                     it_max,
                     fc_max,
                     ecoeff,
                     ccoeff)} {}

CoordSearch::Pattern::Pattern(Index n, Index ptype)
    : pat{vec<std::pair<D, unity_coord>>(
          (check_ptype_and_return(ptype) ==
           static_cast<Index>(Type::GPS2N))
              ? 2 * n
              : n + 1,
          std::pair<D, unity_coord>{
              get_infinity(),
              unity_coord(std::make_pair(vec<Index>(1, 0),
                                         std::vector<D>(1, 0.0)))})},
      poll_type{ptype}
{

  const bool ptype_2n{ptype ==
                      static_cast<Index>(Type::GPS2N)};
  const bool ptype_np1{ptype ==
                       static_cast<Index>(Type::GPSNp1)};

  for (Index i{0}; i < n; ++i) // row
  {
    pat[i].second.first[0] = i;

    if (ptype_2n)
    {
      pat[i + n].second.first[0] = i;
    }

    for (Index j{0}; j < n; ++j) // column
    {
      if (i == j)
      {
        pat[i].second.second[0] = 1.0;
        if (ptype == 0 || ptype == 2)
        {
          pat[i + n].second.second[0] = -1.0;
        }
        break;
      }
    }
  }

  if (ptype_np1)
  {
    pat[n].second.first = vec<Index>(n);
    pat[n].second.second = vec<D>(n);
    for (Index i{0}; i < n; ++i)
    {
      pat[n].second.first[i] = i;
      pat[n].second.second[i] = -1.0;
    }
  }
}

// void CoordSearch::set_pattern(Index n, Index ptype)
// {
//  pattern = Pattern(n, ptype);
// }

bool CoordSearch::Pattern::
    Pattern::less_than::
    operator()(const std::pair<D, unity_coord> &a,
               const std::pair<D, unity_coord> &b)
{
  return a.first < b.first;
}

void CoordSearch::
    new_xiter(const std::vector<D> &xbefore,
              Index i_col,
              std::vector<D> &xnew)
{
  /* computes a new point iterate given step length slen 
  and unit coordinate vector at column i_col */

  for (Index j{0}; j < get_n(); j++) // column
  {
    for (Index i{0}; i < get_row_count(i_col); ++i) // row
    {
      if (pattern.pat[i_col].second.first[i] == j)
      {
        xnew[j] = xbefore[j] + sl * pattern.pat[i_col].second.second[i];
      }
      else
      {
        xnew[j] = xbefore[j];
      }
    }
  }
}

std::ostream &operator<<(std::ostream &os,
                         const CoordSearch::Pattern &p)
{

  Index nn{p.get_col_count()};
  bool np1_ptype{p.poll_type ==
                 static_cast<Index>(
                     CoordSearch::
                         Pattern::
                             Type::GPSNp1)};

  os << std::setw(2) << nn << std::endl;

  if (np1_ptype)
  {
    nn = nn - 1;
  }

  // print  elements of unit direction vector
  for (Index i{0}; i < nn; ++i)
  {
    os << std::setw(2) << p.pat[i].second.second[0];
    if (i < nn - 1)
    {
      os << " ";
    }
  }

  if (np1_ptype)
  {
    os << " [";
    for (Index i{0}; i < p.get_row_count(nn); ++i)
    {
      os << std::setw(2) << p.pat[nn].second.second[i];
      if (i < p.get_row_count(nn) - 1)
      {
        os << ", ";
      }
    }
    os << "]" << std::endl;
  }
  else
  {
    os << std::endl;
  }

  /* print row id */
  for (Index i{0}; i < nn; ++i)
  {
    os << std::setw(2) << p.pat[i].second.first[0];
    if (i < nn - 1)
    {
      os << " ";
    }
  }

  if (np1_ptype)
  {
    os << " [";
    for (Index i{0}; i < p.get_row_count(nn); ++i)
    {
      os << std::setw(2) << p.pat[nn].second.first[i];
      if (i < p.get_row_count(nn) - 1)
      {
        os << ", ";
      }
    }
    os << "]" << std::endl;
  }
  os << std::endl;

  /* print function value */
  if (np1_ptype)
  {
    nn = nn + 1;
  }
  for (Index i{0}; i < nn; ++i)
  {
    os << std::setw(2) << p.pat[i].first;
    if (i < nn - 1)
    {
      os << " ";
    }
  }
  os << std::endl;
  return os;
}

bool CoordSearch::dont_stop()
{
  /* stopping criteria 
  */
  if (sl < par.smin)
  {
    xt_f = 0;
    return false;
  }
  else if (fc >= par.fc_max) // allows for parallel case
  {
    xt_f = 1;
    return false;
  }
  else if (it == par.it_max)
  {
    xt_f = 2;
    return false;
  }
  return true;
}

bool CoordSearch::stop()
{
  return !dont_stop();
}

std::vector<D> CoordSearch::
    copy(const D *x, Index n) const
{
  std::vector<D> v(n);
  std::copy(&x[0], &x[0] + n, std::begin(v));
  return v;
}

void CoordSearch::copy(const D *x,
                       Index n,
                       std::vector<D> &dest)
{
  std::copy(&x[0], &x[0] + n, std::begin(dest));
}

void CoordSearch::
    copy(const std::vector<D> &src,
         std::vector<D> &dest)
{
  std::copy(std::begin(src), std::end(src),
            std::begin(dest));
}

Index CoordSearch::
    check_range_and_get_idx(Index i,
                            Index low,
                            Index high)
{
  if (i < low || i > high)
  {
    Exception::error(
        __FILE__,
        __LINE__,
        std::string{"range error"} +
            std::string{" i = "} +
            std::to_string(i) +
            std::string{" where (low, high) = ("} +
            std::to_string(low) +
            std::string{", "} +
            std::to_string(high) + std::string{")."});
  }
  return i;
}

Index CoordSearch::Pattern::
    check_ptype_and_return(Index ptype)
{
  if (ptype != static_cast<Index>(Type::GPS2N) &&
      ptype != static_cast<Index>(Type::GPSNp1))
  {
    Exception::error(
        __FILE__, __LINE__,
        std::string{"use either 0 (for GPS2N),}"} +
            std::string{" or 1 (for GPSNp1)"});
  }

  return ptype;
}

void CoordSearch::search(std::ostream &os)
{
  print_progress(os);
  init_fval();
  print_progress(os);
  while (!stop())
  {
    poll();
    update_pattern();
    inc_iters();
    print_progress(os);
  }
  print_result(os);
}

bool CoordSearch::poll(Index i, Index &fcalls)
{
  /* poll nearby point at direction i */
  fcalls = 0; /* initialize function evaluation counter */
  if (fc >= par.fc_max || i >= pattern_size())
  {
    return true;
  }

  new_xiter(pt, i, ptnew);

  fxnew = fval(ptnew);
  fcalls++;

  pattern.pat[i].first = fxnew;

  return false;
}

bool CoordSearch::accept_ptnew()
{
  if (fxnew < fxbefore)
  {
    update_fxptbefore();
    return true;
  }
  else
  {
    return false;
  }
}

void CoordSearch::poll()
{
  /* poll nearby points*/

  set_fxptbefore();
  Index i{0};
  Index fcalls{0}; /* function evaluation counter */
  while (!poll(i++, fcalls))
  {
    inc_fcalls(fcalls);

    // if (fxnew < fxbefore)
    // {
    //   update_fxptbefore();
    //   break;
    // }

    if (accept_ptnew())
    {
      break;
    }
  }
  update();
}

void CoordSearch::update_slen()
{
  /* accepts iteration k as successful or not 
  if succesful, updates pattern
  */
  if (fxnew < fx)
  {
    // if successful, accept new point and expand steplength
    fx = fxnew;
    pt = ptnew;
    sl *= par.eco;
  }
  else // if unsuccessful, reduce step length
  {
    sl *= par.cco;
  }
}

void CoordSearch::set_fxptbefore()
{
  fxbefore = fx;
  ptbefore = pt;
}

void CoordSearch::update_fxptbefore()
{
  fxbefore = fxnew;
  ptbefore = ptnew;
}

void CoordSearch::update_fxptnew()
{
  /* update fxnew, and ptnew */
  fxnew = fxbefore;
  ptnew = ptbefore;
}

void CoordSearch::update()
{
  /* update fxnew, and ptnew */
  update_fxptnew();

  /* restore values before */
  set_fxptbefore();

  /* accept new point and update steplength */
  update_slen();
}

void CoordSearch::print_progress(std::ostream &os)
{
  os << std::setw(3) << it << " " << fc << " "
     << sl << " ["
     << std::scientific
     << std::setprecision(6)
     << fx
     << " at (";
  for (size_t j{0}; j < pt.size(); j++)
  {
    os << std::setprecision(6)
       << std::scientific
       << pt[j];
    if (j < pt.size() - 1)
    {
      os << ", ";
    }
    else if (j == pt.size() - 1)
    {
      os << ")]" << std::endl;
    }
  }
}

void CoordSearch::print_result(std::ostream &os)
{
  os << std::endl
     << std::endl
     << std::endl
     << "COORDINATE SEARCH USED "
     << it
     << " ITERATIONS, AND "
     << fc
     << " FUNEVALS, AND RETURNED f(x) = "
     << std::setprecision(8)
     << fx
     << " AT "
     << std::endl;

  for (size_t j{0}; j < pt.size(); j++)
    os << "x[" << std::setw(2)
       << j
       << "] = "
       << std::setprecision(8)
       << std::scientific
       << std::fixed
       << pt[j]
       << std::endl;
  os << "EXIT FLAG: "
     << std::setprecision(8)
     << std::scientific
     << std::fixed
     << xt_f << std::endl
     << "LAST STEPLENGTH: " << sl << std::endl;
}

} // namespace DTDP