/* inc/testfunc.hpp
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

#ifndef test_simplex_hpp
#define test_simplex_hpp

#include <iostream>
#include <cmath>

inline double Rosebrock(const double *x, int n)
{
  /*
    Rosenbrock's parabolic valley
    */
  if (n != 2)
  {
    std::cerr << "In (" << __FILE__
              << " : " << __LINE__
              << ") : number of variables must be 2 but "
              << n << " provided." << std::endl;
    exit(1);
  }

  double fx1 = 100.0 * std::pow(x[1] - x[0] * x[0], 2.0);
  double fx2 = std::pow(1.0 - x[0], 2.0);
  return fx1 + fx2;
}

inline double Powell(const double *x, int n)
{
  /*Powell's quartic function */
  if (n != 4)
  {
    std::cerr << "In (" << __FILE__
              << " : " << __LINE__
              << ") : number of variables must be 4 but "
              << n << " provided." << std::endl;
    exit(1);
  }

  double fx1 = std::pow(x[0] + 10.0 * x[1], 2.0);
  double fx2 = 5.0 * std::pow(x[2] - x[3], 2.0);
  double fx3 = std::pow(x[1] - 2.0 * x[2], 4.0);
  double fx4 = 10.0 * std::pow(x[0] - x[3], 4.0);

  return fx1 + fx2 + fx3 + fx4;
}

inline double Paraboloid(const double *v, int n)
{
  /* Paraboloid center at (1,2), scale factors (10, 20),
   minimum value 30 */

  if (n != 2)
  {
    std::cerr << "In (" << __FILE__
              << " : " << __LINE__
              << ") : number of variables must be 2 but "
              << n << " provided." << std::endl;
    exit(1);
  }

  double p[5] = {1.0, 2.0, 10.0, 20.0, 30.0};

  double x = v[0];
  double y = v[1];

  return p[2] * (x - p[0]) * (x - p[0]) +
         p[3] * (y - p[1]) * (y - p[1]) + p[4];
}

#endif