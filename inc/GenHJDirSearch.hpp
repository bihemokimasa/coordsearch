/* inc/GenHJDirSearch.hpp
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

//
// I extend the baseline Hooke and Jeeves algorithm
// to solve the general constrained optimization problem
//
// min f(x) s.t. h(x) <= 0
//
// with bound constraints such that LB <= x <= UB
//
// where h(x) is the set of all relaxable constraints

#ifndef GENHJDIRSEARCH_HPP
#define GENHJDIRSEARCH_HPP

#include "HJDirSearch.hpp"
#include "Array1D.hpp"
#include "SolverParam.hpp"
// #include "HJParams.hpp"


namespace DTDP
{
class GenHJDirSearch : public HJDirSearch
{
public:
  using D = DTDP::D;
  using Index = DTDP::Index;
  // const D TOL{1E-6};

  GenHJDirSearch() = default;
  GenHJDirSearch(const GenHJDirSearch &) = default;
  GenHJDirSearch &operator=(const GenHJDirSearch &) = default;
  GenHJDirSearch(GenHJDirSearch &&) = default;
  GenHJDirSearch &operator=(GenHJDirSearch &&) = default;
  ~GenHJDirSearch() = default;

  // H is a function pointer to the set of constraints functions
  // that returns an array of size equal to the
  // dimensions of aIneqCons
  GenHJDirSearch(const HJDirSearch &hj,
                 const std::function<void(
                     const Array1D<D> &x,
                     Array1D<D> &constr_result,
                     const SolverParam &params)> &hh,
                 const Index countIneqCons,
                 const D tol = 1E-6);

  GenHJDirSearch(const HJDirSearch &hj, const Array1D<D> &lowerbound,
                 const Array1D<D> &upperbound, const D tol = 1E-6);

  GenHJDirSearch(const HJDirSearch &hj,
                 const std::function<void(const Array1D<D> &x,
                                          Array1D<D> &constr_result,
                                          const SolverParam &params)> &hh,
                 const Array1D<D> &lowerbound,
                 const Array1D<D> &upperbound,
                 const Index countIneqCons = 0,
                 const D tol = 1E-6);

  GenHJDirSearch(HJDirSearch &&hj,
                 const std::function<void(
                     const Array1D<D> &x,
                     Array1D<D> &constr_result, const SolverParam &params)> &hh,
                 const Index countIneqCons,
                 const D tol = 1E-6);

  GenHJDirSearch(HJDirSearch &&hj,
                 const std::function<void(
                     const Array1D<D> &x,
                     Array1D<D> &constr_result, const SolverParam &params)> &hh,
                 const Array1D<D> &lowerbound,
                 const Array1D<D> &upperbound,
                 const Index countIneqCons = 0,
                 const D tol = 1E-6);

  // constrained violation function:
  // returns the sum of the square of the violations of each relaxable constraints
  D MaxConsViol(const Array1D<D> &point, const SolverParam &params);
  void SetMaxConsViol(const Array1D<D> &point, const SolverParam &params)
  {
    aMaxConsViol = MaxConsViol(point, params);
  }
  void SetMaxConsViol(D val) { aMaxConsViol = val; }
  void AdjustToBounds(Array1D<D> &point);

  void SetLowerBound(const Array1D<D> &lb) { aLowerBound = lb; }
  void SetLowerBound(const D lb, const Index i) { aLowerBound(i) = lb; }
  const Array1D<D> &LowerBound() const { return aLowerBound; }

  void SetUpperBound(const Array1D<D> &ub) { aUpperBound = ub; }
  void SetUpperBound(const D ub, const Index i) { aUpperBound(i) = ub; }
  const Array1D<D> &UpperBound() const { return aUpperBound; }

  void SetIneqCons(const Array1D<D> &con)
  {
    aIneqCons = con;
    aConstArr.constArr = Array1D<D>(con.size());
  }
  const Array1D<D> &IneqCons() const { return aIneqCons; }

  void SetTol(const D tol) { aTol = tol; }
  D Tol() const { return aTol; }

  D MaxConsViol() const { return aMaxConsViol; }

  void SetConsFun(const std::function<void(
                      const Array1D<D> &x,
                      Array1D<D> &constr_result,
                      const SolverParam &params)> &cf)
  {
    aConsFun = cf;
  }

  bool SatisfiesCons(const Array1D<D> &pt,
                     const SolverParam &params)
  {
    return MaxConsViol(pt, params) <= aTol;
  }

  // Evaluate a point and where constraints are violated return infinity
  // this is the extreme barrier method
  virtual D EvalPt(Array1D<D> &pt,
                   const SolverParam &params) override;
  virtual int hooke(const SolverParam &params,
                    std::ostream &os = std::cout) override;
  virtual double best_nearby(Array1D<D> &delta,
                             Array1D<D> &point,
                             int &funevals,
                             D prevbest,
                             const SolverParam &params) override;

  // the result of each of the constraints in this array
  struct ConstArr
  {
    Array1D<D> constArr{Array1D<D>()};
  } aConstArr;

private:
  Array1D<D> aLowerBound{Array1D<D>()};
  Array1D<D> aUpperBound{Array1D<D>()};
  Array1D<D> aIneqCons{Array1D<D>()}; // value of the inequality constraint
  D aTol{1E-6};                       // tolenrance on constraint violation
  D aMaxConsViol{0.0};
  std::function<void(const Array1D<D> &x,
                     Array1D<D> &constr_result,
                     const SolverParam &params)>
      aConsFun{nullptr};
};
} // namespace DTDP

#endif