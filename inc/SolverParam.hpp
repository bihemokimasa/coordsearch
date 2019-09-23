/* inc/SolverParams.hpp
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
* Can pass parameters to the solver as instance of struct SolverParam
*/

#ifndef SOLVERPARAM_HPP
#define SOLVERPARAM_HPP

#include "Array1D.hpp"

namespace DTDP
{
struct SolverParam
{
  D p; 
};
} // namespace DTDP

#endif