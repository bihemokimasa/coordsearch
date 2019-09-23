/* src/Exception.cpp
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

#include "Exception.hpp"

#include <errno.h>
#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>
#include <iomanip>

namespace DTDP
{

Exception::Exception(std::string file,
                     int line,
                     std::string msg,
                     int error_number)
    : std::runtime_error(std::move(msg)),
      aFile{std::move(file)},
      aLine{line},
      aErrorNumber{error_number} {}

const char *Exception::what(void) const throw()
{
  std::ostringstream oss;
  if (!aFile.empty() || aLine > 0)
    oss << "DTDP::Exception thrown"
        << " (" << aFile << ", " << aLine << ")"
        << ": " << std::runtime_error::what()
        << std::endl;

  aWhat = oss.str();
  return aWhat.c_str();
}

void Exception::error(std::string file,
                      int line,
                      std::string msg,
                      int errnum)
{
  throw(Exception(std::move(file),
                  line,
                  std::move(msg),
                  errnum));
}

// Courtesy of: http://man7.org/linux/man-pages/man3/backtrace.3.html
void Exception::PrintStackTrace()
{
  int j, nptrs;
  int BT_BUF_SIZE{100};
  void *buffer[BT_BUF_SIZE];
  char **strings;

  nptrs = backtrace(buffer, BT_BUF_SIZE);
  printf("backtrace() returned %d addresses\n", nptrs);

  /* The call backtrace_symbols_fd(buffer, nptrs, STDOUT_FILENO)
              would produce similar output to the following: */

  strings = backtrace_symbols(buffer, nptrs);
  if (strings == NULL)
  {
    perror("backtrace_symbols");
    exit(EXIT_FAILURE);
  }

  for (j = 0; j < nptrs; j++)
  {
    printf("[bt, %2d]: %s\n", j, strings[j]);

    // Courtesy of: https://stackoverflow.com/questions/3151779/best-way-to-invoke-gdb-from-inside-program-to-print-its-stacktrace/4611112#4611112
    /* find first occurence of '(' or ' ' in message[i] and assume
     * everything before that is the file name. (Don't go beyond 0 though
     * (string terminator)*/
    int p = 0;
    while (strings[j][p] != '(' && strings[j][p] != ' ' && strings[j][p] != 0)
      ++p;

    char syscom[256];
    sprintf(syscom, "addr2line %p -e %.*s", buffer[j], p, strings[j]);
    //last parameter is the file name of the symbol
    system(syscom);
  }

  free(strings);
}
} // namespace DTDP
