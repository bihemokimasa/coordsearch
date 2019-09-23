/* inc/Exception.hpp
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
* A custom Exception class, for exception handling
*/

#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>

namespace DTDP
{

enum class ERROR_NUMBER
{
  IGNORE = 0,  // error number not provided, ingore
  TOO_FEW_OBS, // too few observations
  NAN_TO_ZERO, // when we set NaN value to zero
  NONREAL,     // non-real result
  RANGE_ERROR,
};

class Exception : virtual public std::runtime_error
{
public:
  Exception() = default;
  Exception(const Exception &) = default;
  Exception &operator=(const Exception &) = default;

  explicit Exception(std::string file,
                     int line,
                     std::string msg,
                     int error_number);

  virtual ~Exception(void) throw() {}

  const char *what(void) const throw() override;

  int ErrorNumber() const { return aErrorNumber; }

  static void error(std::string file,
                    int line,
                    std::string msg,
                    int error_number = static_cast<int>(ERROR_NUMBER::IGNORE));

  static void PrintStackTrace();

private:
  mutable std::string aWhat{std::string{}}; // Error message.
  std::string aFile{std::string{}};         // File where the exception is thrown.
  int aLine{0};                             // Line number at which the exception is thrown.
  int aErrorNumber{0};
  // bool aPrintStackTrace{false}; // if true prints backtrace
};

} // namespace DTDP

#endif