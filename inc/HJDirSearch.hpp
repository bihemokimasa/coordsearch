/* inc/HJDirSearch.hpp
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

// a Wrapper for the Hooke direct search algorithm

#ifndef HJDIRSEARCH_HPP
#define HJDIRSEARCH_HPP



#include <functional>

#include "SolverParam.hpp"
#include "Array1D.hpp"

namespace DTDP
{
class HJDirSearch
{
public:
  using D = DTDP::D;
  using Index = DTDP::Index;
  static constexpr D default_rho{0.5};
  static constexpr D default_epsilon{1E-6};

  HJDirSearch() = default;
  HJDirSearch(const HJDirSearch &) = default;
  HJDirSearch &operator=(const HJDirSearch &) = default;
  HJDirSearch(HJDirSearch &&) = default;
  HJDirSearch &operator=(HJDirSearch &&) = default;
  virtual ~HJDirSearch() = default;

  HJDirSearch(const Array1D<D> &startpt,
              const std::function<double(
                  const Array1D<D> &x,
                  const SolverParam &params)> &ff);

  HJDirSearch(const Array1D<D> &startpt,
              const std::function<double(
                  const Array1D<D> &x,
                  const SolverParam &params)> &ff,
              D rho,
              D eps,
              Index itermax);

  virtual int hooke(const SolverParam &params,
                    std::ostream &os = std::cout);
  virtual D EvalPt(Array1D<D> &pt,
                   const SolverParam &params)
  {
    return aObjFun(pt, params);
  } // evaluate point


  virtual void PrintEndMsg(std::ostream &os) const;

  void SetEndPt(const Array1D<D> &pt) { aEndPt = pt; }
  void SetEndPt(const D pt, const Index i) { aEndPt(i) = pt; }
  const Array1D<D> &EndPt() const { return aEndPt; }
  D EndPt(const Index i) const { return aEndPt(i); }

  void SetRho(D rho = default_rho) { aRho = rho; }
  D Rho() const { return aRho; }

  void SetEpsilon(D eps = default_epsilon) { aEpsilon = eps; }
  D Epsilon() const { return aEpsilon; }

  D StepLength() const { return aStepLength; }
  void SetStepLength(D val) { aStepLength = val; }
  D MinVal() const { return aMinVal; }
  void SetMinVal(D val) { aMinVal = val; }

  Index FunEvals() const { return aFunEvals; }
  void SetFunEvals(Index fval) { aFunEvals = fval; }
  void SetMaxFunEvals(Index val) { aMaxFunEvals = val; }

  Index DefaultMaxFunEvals() const
  {
    return (aEndPt.size() == 0)
               ? 1000000
               : aEndPt.size() * 2000;
  };

  Index MaxFunEvals() const { return aMaxFunEvals; }

  Index Iters() const { return aIters; }
  void SetIters(Index val) { aIters = val; }
  Index IterMax() const { return aIterMax; }

  Index DefaultIterMax() const
  {
    return (aEndPt.size() == 0)
               ? 1000
               : aEndPt.size() * 100;
  };

  void SetIterMax(const Index maxiter) { aIterMax = maxiter; }
  void SetDisplay(const bool &disp) { aDisplay = disp; }
  bool Display() const { return aDisplay; }

  Index ExitFlag() const { return aExitFlag; }
  void SetExitFlag(Index val) { aExitFlag = val; }

  void SetObjFun(const std::function<double(
                     const Array1D<D> &x,
                     const SolverParam &params)> &ff)
  {
    aObjFun = ff;
  }

  virtual bool Stop(D steplength,
                    Index iters,
                    Index funevals)
  {
    if (steplength <= aEpsilon)
    {
      aExitFlag = 0;
      return true;
    }
    else if ((iters >= aIterMax) ||
             (funevals >= aMaxFunEvals))
    {
      aExitFlag = 1;
      return true;
    }
    return false;
  }

  int hooke(double &fbefore,
            const Array1D<D> &startpt,
            Array1D<D> &endpt,
            double &steplength,
            int &iters,
            int &funevals,
            D rho,
            D epsilon,
            const SolverParam &params,
            std::ostream &os,
            const bool &display = true);

  virtual double best_nearby(Array1D<D> &delta,
                             Array1D<D> &point, int &funevals,
                             D prevbest,
                             const SolverParam &params);

  virtual void PrintProgress(Index funevals, D fbefore, const Array1D<D> &xbefore, std::ostream &os) const;


  // const D RHO{0.5};          // stepsize geometric shrink
  // const D EPS{1E-6};         //  ending value of stepsize
  // const Index MAXITER{5000}; // max # of iterations
private:
  Array1D<D> aEndPt{Array1D<D>()}; //  This is the location of  the local minimum, calculated by the program
  D aRho{0.5};                     // the algorithm convergence control
  D aEpsilon{1E-5};                // This is the criterion for halting the search for a minimum
  D aStepLength{0.0};
  D aMinVal{0.0};
  Index aFunEvals{0};
  Index aMaxFunEvals{0};
  Index aIters{0};
  Index aIterMax{5000}; // A second, rarely used, halting  criterion.  If the algorithm uses >= itermax iterations, halt
  bool aDisplay{false};
  Index aExitFlag{0};
  std::function<double(
      const Array1D<D> &x,
      const SolverParam &params)>
      aObjFun{nullptr};
  // D(*ObjFun)
  // (const Array1D<D> &x){nullptr};
};
} // namespace DTDP
#endif