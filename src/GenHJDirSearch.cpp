/* src/GenHJDirSearch.cpp
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

#include "GenHJDirSearch.hpp"

#include <limits>
#include <cassert>

namespace DTDP
{
  GenHJDirSearch::GenHJDirSearch(const HJDirSearch &hj,
                                 const std::function<void(
                                     const Array1D<D> &x,
                                     Array1D<D> &constr_result,
                                     const SolverParam &params)> &hh,
                                 const Index countIneqCons,
                                 const D tol)
      : HJDirSearch{hj},
        aLowerBound{Array1D<D>(hj.EndPt().size())},
        aUpperBound{{Array1D<D>(hj.EndPt().size())}},
        aIneqCons{Array1D<D>(countIneqCons)},
        aTol{tol},
        aConsFun{hh}
  {
    for (Index i{0}; i < hj.EndPt().size(); i++)
    {
      aLowerBound(i) = -std::numeric_limits<double>::infinity();
      aUpperBound(i) = std::numeric_limits<double>::infinity();
    }
    aConstArr.constArr = Array1D<D>(countIneqCons);
  }

  GenHJDirSearch::GenHJDirSearch(const HJDirSearch &hj,
                                 const Array1D<D> &llb,
                                 const Array1D<D> &uub, const D tol)
      : HJDirSearch{hj},
        aLowerBound{llb},
        aUpperBound{uub},
        aTol{tol} {}

  GenHJDirSearch::GenHJDirSearch(const HJDirSearch &hj,
                                 const std::function<void(
                                     const Array1D<D> &x,
                                     Array1D<D> &constr_result,
                                     const SolverParam &params)> &hh,
                                 const Array1D<D> &llb,
                                 const Array1D<D> &uub,
                                 const Index countIneqCons,
                                 const D tol)
      : HJDirSearch{hj},
        aLowerBound{llb},
        aUpperBound{uub},
        aIneqCons{Array1D<D>(countIneqCons)},
        aTol{tol},
        aConsFun{hh}
  {
    assert((llb.size() == uub.size()) &&
           (llb.size() == hj.EndPt().size()));
    aConstArr.constArr = Array1D<D>(countIneqCons);
  }

  GenHJDirSearch::GenHJDirSearch(HJDirSearch &&hj,
                                 const std::function<void(
                                     const Array1D<D> &x,
                                     Array1D<D> &constr_result,
                                     const SolverParam &params)> &hh,
                                 const Index countIneqCons,
                                 const D tol)
      : HJDirSearch{std::move(hj)},
        aLowerBound{Array1D<D>(this->EndPt().size())},
        aUpperBound{{Array1D<D>(this->EndPt().size())}},
        aIneqCons{Array1D<D>(countIneqCons)},
        aTol{tol},
        aConsFun{hh}
  {
    for (Index i{0}; i < EndPt().size(); i++)
    {
      aLowerBound(i) = -std::numeric_limits<double>::infinity();
      aUpperBound(i) = std::numeric_limits<double>::infinity();
    }
    aConstArr.constArr = Array1D<D>(countIneqCons);
  }

  GenHJDirSearch::GenHJDirSearch(HJDirSearch &&hj,
                                 const std::function<void(
                                     const Array1D<D> &x,
                                     Array1D<D> &constr_result,
                                     const SolverParam &params)> &hh,
                                 const Array1D<D> &llb,
                                 const Array1D<D> &uub,
                                 const Index countIneqCons,
                                 const D tol)
      : HJDirSearch{std::move(hj)},
        aLowerBound{llb},
        aUpperBound{uub},
        aIneqCons{Array1D<D>(countIneqCons)},
        aTol{tol},
        aConsFun{hh}
  {
    assert((llb.size() == uub.size()) &&
           (llb.size() == EndPt().size()));
    aConstArr.constArr = Array1D<D>(countIneqCons);
  }

  D GenHJDirSearch::MaxConsViol(const Array1D<D> &point,
                                const SolverParam &params)
  {
    aConsFun(point, aConstArr.constArr, params);
    D maxConsViol = 0.0;
    // for (Index i{0}; i < aConstArr.constArr.size(); i++)
    //   maxConsViol += std::pow(std::max(aConstArr.constArr(i), 0.0), 2.0);
    // aMaxConsViol += std::max(aConstArr.constArr(i), 0.0);
    // std:: cout << "aMaxConsViol" << aMaxConsViol << std::endl;

    for (Index i{0}; i < aConstArr.constArr.size(); i++)
      maxConsViol = std::max(aConstArr.constArr(i), maxConsViol);
    return maxConsViol;
  }

  void GenHJDirSearch::AdjustToBounds(Array1D<D> &point)
  {
    for (Index i{0}; i < aLowerBound.size(); i++)
    {
      if (point(i) < aLowerBound(i))
        point(i) = aLowerBound(i);
      else if (point(i) > aUpperBound(i))
        point(i) = aUpperBound(i);
    }
  }

  D GenHJDirSearch::EvalPt(Array1D<D> &point,
                           const SolverParam &params)
  {
    if (aConsFun == nullptr)
      return HJDirSearch::EvalPt(point, params);
    if (!SatisfiesCons(point, params))
    {
      SetFunEvals(FunEvals() - 1); // decrement funevals, so that actual funevals are used
      return std::numeric_limits<double>::infinity();
    }
    else
      return HJDirSearch::EvalPt(point, params);
  }

  int GenHJDirSearch::hooke(const SolverParam &params,
                            std::ostream &os)
  {
    if (aConsFun != nullptr)
    {
      aConsFun(EndPt(), aConstArr.constArr, params);
      for (Index i{0}; i < aConstArr.constArr.size(); i++)
        if (aConstArr.constArr(i) > aTol)
        {
          // std::cout << EndPt() << std::endl;
          Exception::error(__FILE__, __LINE__,
                           std::string{"Initial supplied point "} +
                               std::string{"does not satisfy constraints."},
                           static_cast<int>(ERROR_NUMBER::IGNORE));
        }
    }
    Index iters{HJDirSearch::hooke(params, os)};
    if (aConsFun != nullptr)
    {
      SetMaxConsViol(EndPt(), params);
    }
    return iters;
  }

  double GenHJDirSearch::best_nearby(
      Array1D<D> &delta, Array1D<D> &point, int &funevals,
      D prevbest, const SolverParam &params)
  // 		double delta[VARS],
  // 		point[VARS];
  // double prevbest;
  // int nvars;
  {
    // double z[VARS];
    int nvars{point.size()};
    Array1D<D> z{Array1D<D>(nvars)};
    double minf, ftmp;
    int i;
    minf = prevbest;
    AdjustToBounds(point);

    if (funevals < MaxFunEvals())
    {
      for (i = 0; i < nvars; i++)
        z[i] = point[i];

      for (i = 0; i < nvars; i++)
      {
        z[i] = point[i] + delta[i];
        AdjustToBounds(z);
        ftmp = EvalPt(z, params);
        funevals++;
        if (ftmp < minf)
          minf = ftmp;
        else
        {
          delta[i] = 0.0 - delta[i];
          z[i] = point[i] + delta[i];
          AdjustToBounds(z);
          ftmp = EvalPt(z, params);
          funevals++;
          if (ftmp < minf)
            minf = ftmp;
          else
            z[i] = point[i];
        }
      }
      for (i = 0; i < nvars; i++)
        point[i] = z[i];
    }

    return (minf);
  }
} // namespace DTDP
