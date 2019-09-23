/* src/HJDirSearch.cpp
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

#include "HJDirSearch.hpp"

#include <iostream>
#include <cmath>

#include <mpi.h>

#include "Exception.hpp"

namespace DTDP
{
HJDirSearch::
    HJDirSearch(const Array1D<D> &startpt,
                const std::function<double(
                    const Array1D<D> &x,
                    const SolverParam &params)> &ff)
    : aEndPt{startpt},
      aRho{0.5},
      aEpsilon{1E-6},
      aStepLength{0.0},
      aMinVal{0.0},
      aFunEvals{0},
      aMaxFunEvals{DefaultMaxFunEvals()},
      aIters{0},
      aIterMax{DefaultIterMax()},
      aDisplay{false},
      aExitFlag{0},
      aObjFun{ff}
{
}

HJDirSearch::
    HJDirSearch(const Array1D<D> &startpt,
                const std::function<double(
                    const Array1D<D> &x,
                    const SolverParam &params)> &ff,
                D rho,
                D eps,
                Index itermax)
    : HJDirSearch(startpt, ff)
{
  aRho = rho;
  aEpsilon = eps;
  aIterMax = itermax;
}

int HJDirSearch::
    hooke(const SolverParam &params,
          std::ostream &os)
{
  aMinVal = EvalPt(aEndPt, params);

  if (std::isnan(aMinVal) || std::isinf(aMinVal))
    Exception::error(__FILE__, __LINE__,
                     std::string{"Error: Initial function value "} +
                         std::to_string(aMinVal),
                     static_cast<int>(ERROR_NUMBER::IGNORE));

  aFunEvals = 0;
  hooke(aMinVal, aEndPt, aEndPt, aStepLength, aIters,
        aFunEvals, aRho, aEpsilon, params, os, aDisplay);

  if (aDisplay)
    PrintEndMsg(os);
  return aIters;
};

void HJDirSearch::PrintEndMsg(std::ostream &os) const
{
  os << std::endl
     << std::endl
     << std::endl
     << "HOOKE USED "
     << aIters
     << " ITERATIONS, AND "
     << aFunEvals
     << " FUNEVALS, AND RETURNED f(x) = "
     << std::setprecision(8)
     << aMinVal
     << " AT "
     << std::endl;
  for (Index j{0}; j < aEndPt.size(); j++)
    os << "x[" << std::setw(2)
       << j
       << "] = "
       << std::setprecision(4)
       << std::scientific
       << std::fixed
       << aEndPt(j)
       << std::endl;
  os << std::endl;
}

// Includes minor modifications by Bihemo Kimasa

/* Nonlinear Optimization using the algorithm of Hooke and Jeeves  */
/*	12 February 1994	author: Mark G. Johnson 	   */

/* Find a point X where the nonlinear function f(X) has a local    */
/* minimum.  X is an n-vector and f(X) is a scalar.  In mathe-	   */
/* matical notation  f: R^n -> R^1.  The objective function f()    */
/* is not required to be continuous.  Nor does f() need to be	   */
/* differentiable.  The program does not use or require 	   */
/* derivatives of f().						   */

/* The software user supplies three things: a subroutine that	   */
/* computes f(X), an initial "starting guess" of the minimum point */
/* X, and values for the algorithm convergence parameters.  Then   */
/* the program searches for a local minimum, beginning from the    */
/* starting guess, using the Direct Search algorithm of Hooke and  */
/* Jeeves.							   */

/* This C program is adapted from the Algol pseudocode found in    */
/* "Algorithm 178: Direct Search" by Arthur F. Kaupe Jr., Commun-  */
/* ications of the ACM, Vol 6. p.313 (June 1963).  It includes the */
/* improvements suggested by Bell and Pike (CACM v.9, p. 684, Sept */
/* 1966) and those of Tomlin and Smith, "Remark on Algorithm 178"  */
/* (CACM v.12).  The original paper, which I don't recommend as    */
/* highly as the one by A. Kaupe, is:  R. Hooke and T. A. Jeeves,  */
/* "Direct Search Solution of Numerical and Statistical Problems", */
/* Journal of the ACM, Vol. 8, April 1961, pp. 212-229. 	   */

/* Calling sequence:						   */
/*  int hooke(nvars, startpt, endpt, rho, epsilon, itermax)	   */
/*								   */
/*     nvars	   {an integer}  This is the number of dimensions  */
/*		   in the domain of f().  It is the number of	   */
/*		   coordinates of the starting point (and the	   */
/*		   minimum point.)				   */
/*     startpt	   {an array of doubles}  This is the user-	   */
/*		   supplied guess at the minimum.		   */
/*     endpt	   {an array of doubles}  This is the location of  */
/*		   the local minimum, calculated by the program    */
/*     rho	   {a double}  This is a user-supplied convergence */
/*		   parameter (more detail below), which should be  */
/*		   set to a value between 0.0 and 1.0.	Larger	   */
/*		   values of rho give greater probability of	   */
/*		   convergence on highly nonlinear functions, at a */
/*		   cost of more function evaluations.  Smaller	   */
/*		   values of rho reduces the number of evaluations */
/*		   (and the program running time), but increases   */
/*		   the risk of nonconvergence.	See below.	   */
/*     epsilon	   {a double}  This is the criterion for halting   */
/*		   the search for a minimum.  When the algorithm   */
/*		   begins to make less and less progress on each   */
/*		   iteration, it checks the halting criterion: if  */
/*		   the stepsize is below epsilon, terminate the    */
/*		   iteration and return the current best estimate  */
/*		   of the minimum.  Larger values of epsilon (such */
/*		   as 1.0e-4) give quicker running time, but a	   */
/*		   less accurate estimate of the minimum.  Smaller */
/*		   values of epsilon (such as 1.0e-7) give longer  */
/*		   running time, but a more accurate estimate of   */
/*		   the minimum. 				   */
/*     itermax	   {an integer}  A second, rarely used, halting    */
/*		   criterion.  If the algorithm uses >= itermax    */
/*		   iterations, halt.				   */

/* The user-supplied objective function f(x,n) should return a C   */
/* "double".  Its  arguments are  x -- an array of doubles, and    */
/* n -- an integer.  x is the point at which f(x) should be	   */
/* evaluated, and n is the number of coordinates of x.	That is,   */
/* n is the number of coefficients being fitted.		   */

/* rho, the algorithm convergence control			   */
/*	The algorithm works by taking "steps" from one estimate of */
/*    a minimum, to another (hopefully better) estimate.  Taking   */
/*    big steps gets to the minimum more quickly, at the risk of   */
/*    "stepping right over" an excellent point.  The stepsize is   */
/*    controlled by a user supplied parameter called rho.  At each */
/*    iteration, the stepsize is multiplied by rho  (0 < rho < 1), */
/*    so the stepsize is successively reduced.			   */
/*	Small values of rho correspond to big stepsize changes,    */
/*    which make the algorithm run more quickly.  However, there   */
/*    is a chance (especially with highly nonlinear functions)	   */
/*    that these big changes will accidentally overlook a	   */
/*    promising search vector, leading to nonconvergence.	   */
/*	Large values of rho correspond to small stepsize changes,  */
/*    which force the algorithm to carefully examine nearby points */
/*    instead of optimistically forging ahead.	This improves the  */
/*    probability of convergence.				   */
/*	The stepsize is reduced until it is equal to (or smaller   */
/*    than) epsilon.  So the number of iterations performed by	   */
/*    Hooke-Jeeves is determined by rho and epsilon:		   */
/*	    rho**(number_of_iterations) = epsilon		   */
/*	In general it is a good idea to set rho to an aggressively */
/*    small value like 0.5 (hoping for fast convergence).  Then,   */
/*    if the user suspects that the reported minimum is incorrect  */
/*    (or perhaps not accurate enough), the program can be run	   */
/*    again with a larger value of rho such as 0.85, using the	   */
/*    result of the first minimization as the starting guess to    */
/*    begin the second minimization.				   */

/* Normal use: (1) Code your function f() in the C language	   */
/*	       (2) Install your starting guess {or read it in}	   */
/*	       (3) Run the program				   */
/*	       (4) {for the skeptical}: Use the computed minimum   */
/*		      as the starting point for another run	   */

/* Data Fitting:						   */
/*	Code your function f() to be the sum of the squares of the */
/*	errors (differences) between the computed values and the   */
/*	measured values.  Then minimize f() using Hooke-Jeeves.    */
/*	EXAMPLE: you have 20 datapoints (ti, yi) and you want to   */
/*	find A,B,C such that  (A*t*t) + (B*exp(t)) + (C*tan(t))    */
/*	fits the data as closely as possible.  Then f() is just    */
/*	f(x) = SUM (measured_y[i] - ((A*t[i]*t[i]) + (B*exp(t[i])) */
/*				  + (C*tan(t[i]))))^2		   */
/*	where x[] is a 3-vector consisting of {A, B, C}.	   */

/*								   */
/*  The author of this software is M.G. Johnson.		   */
/*  Permission to use, copy, modify, and distribute this software  */
/*  for any purpose without fee is hereby granted, provided that   */
/*  this entire notice is included in all copies of any software   */
/*  which is or includes a copy or modification of this software   */
/*  and in all copies of the supporting documentation for such	   */
/*  software.  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT    */
/*  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE   */
/*  AUTHOR NOR AT&T MAKE ANY REPRESENTATION OR WARRANTY OF ANY	   */
/*  KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS    */
/*  FITNESS FOR ANY PARTICULAR PURPOSE. 			   */
/*								   */

int HJDirSearch::
    hooke(double &fbefore,
          const Array1D<D> &startpt,
          Array1D<D> &endpt,
          double &steplength,
          int &iters,
          int &funevals,
          D rho,
          D epsilon,
          const SolverParam &params,
          std::ostream &os,
          const bool &display)
{
  // double delta[VARS];
  int nvars{startpt.size()};
  Array1D<D> delta{Array1D<D>(nvars)};
  double newf, tmp;
  Array1D<D> xbefore{Array1D<D>(nvars)},
      newx{Array1D<D>(nvars)};
  // double xbefore[VARS], newx[VARS];
  int i, keep;
  int iadj;
  for (i = 0; i < nvars; i++)
  {
    newx[i] = xbefore[i] = startpt[i];
    delta[i] = fabs(startpt[i] * rho);
    if (delta[i] == 0.0)
      delta[i] = rho;
  }
  iadj = 0;
  steplength = rho;
  iters = 0;
  // fbefore = aMinVal;
  funevals++;
  newf = fbefore;
  bool stop = false;
  while (stop == false)
  {
    iters++;
    iadj++;
    if (display)
      PrintProgress(funevals, fbefore, xbefore, os);

    /* find best new point, one coord at a time */
    for (i = 0; i < nvars; i++)
    {
      newx[i] = xbefore[i];
    }
    newf = best_nearby(delta,
                       newx,
                       funevals,
                       fbefore,
                       params);
    /* if we made some improvements, pursue that direction */
    keep = 1;
    while ((newf < fbefore) && (keep == 1))
    {
      iadj = 0;
      for (i = 0; i < nvars; i++)
      {
        /* firstly, arrange the sign of delta[] */
        if (newx[i] <= xbefore[i])
          delta[i] = 0.0 - fabs(delta[i]);
        else
          delta[i] = fabs(delta[i]);
        /* now, move further in this direction */
        tmp = xbefore[i];
        xbefore[i] = newx[i];
        newx[i] = newx[i] + newx[i] - tmp;
      }
      fbefore = newf;
      newf = best_nearby(delta,
                         newx,
                         funevals,
                         fbefore,
                         params);

      /* if the further (optimistic) move was bad.... */
      if (newf >= fbefore)
        break;
      /* make sure that the differences between the new */
      /* and the old points are due to actual */
      /* displacements; beware of roundoff errors that */
      /* might cause newf < fbefore */
      keep = 0;
      for (i = 0; i < nvars; i++)
      {
        keep = 1;
        if (fabs(newx[i] - xbefore[i]) >
            (0.5 * fabs(delta[i])))
          break;
        else
          keep = 0;
      }
    }
    if ((steplength >= epsilon) && (newf >= fbefore))
    {
      steplength = steplength * rho;
      for (i = 0; i < nvars; i++)
      {
        delta[i] *= rho;
      }
    }

    // (iters < itermax) &&
    //      (steplength > epsilon) &&
    //      (funevals < aMaxFunEvals)
    stop = Stop(steplength, iters, funevals);
  }
  for (i = 0; i < nvars; i++)
    endpt[i] = xbefore[i];
  return (iters);
}

double HJDirSearch::best_nearby(
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

  if (funevals < MaxFunEvals())
  {
    for (i = 0; i < nvars; i++)
      z[i] = point[i];
    for (i = 0; i < nvars; i++)
    {
      z[i] = point[i] + delta[i];
      ftmp = EvalPt(z, params);
      funevals++;
      if (ftmp < minf)
        minf = ftmp;
      else
      {
        delta[i] = 0.0 - delta[i];
        z[i] = point[i] + delta[i];
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

void HJDirSearch::PrintProgress(Index funevals,
                                D fbefore,
                                const Array1D<D> &xbefore,
                                std::ostream &os) const
{
  os << std::endl
     << "After " << std::setw(5)
     << funevals
     << " funevals, f(x) = "
     << std::scientific
     << std::fixed
     << std::setprecision(4)
     << fbefore
     << " at "
     << std::endl;
  for (Index j = 0; j < xbefore.size(); j++)
    os << "   x["
       << std::setw(2)
       << j
       << "] = "
       << std::scientific
       << std::fixed
       << std::setprecision(4)
       << xbefore[j]
       << std::endl;
}
} // namespace DTDP