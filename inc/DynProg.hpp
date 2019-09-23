/* inc/DynProg.hpp
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
* Provides some structs and consts for use throughout the programm
* For instance, these structs are the convergence parameters or names of files
*/ 



#ifndef DYNGRID_H
#define DYNGRID_H

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>

#include <mpi.h>

#define MASTER_RANK 0
#define DOUBLE_PRECISION 16

namespace DTDP
{

using Index = int;
using D = double;
const double PI = std::acos(-1);
const double Inf = std::numeric_limits<double>::infinity();

} // namespace DTDP

#endif // DynGrid.hpp
