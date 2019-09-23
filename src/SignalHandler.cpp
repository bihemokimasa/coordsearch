/* src/SignalHandler.cpp
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

#include "SignalHandler.hpp"

#include <signal.h>
#include <errno.h>
#include <execinfo.h>
#include <stdlib.h>
#include <unistd.h>
#include <ucontext.h>

#include "Exception.hpp"

namespace DTDP
{

bool SignalHandler::aGotExitSignal = false;

bool SignalHandler::GotExitSignal()
{
  return aGotExitSignal;
}

void SignalHandler::SetExitSignal(bool exitsignal)
{
  aGotExitSignal = exitsignal;
}

void SignalHandler::ExitSignalHandler(int signum)
{
  aGotExitSignal = true;
  SetExitSignal(true);

  // Courtesy of http://www.adv-ci.com/blog/2012/12/27/cc-stack-trace/
  // associate each signal with a signal name string.
  const char *name = NULL;

  // get signal name
  switch (signum)
  {
  case SIGABRT:
    name = "SIGABRT";
    break;
  case SIGSEGV:
    name = "SIGSEGV";
    break;
  case SIGBUS:
    name = "SIGBUS";
    break;
  case SIGILL:
    name = "SIGILL";
    break;
  case SIGFPE:
    name = "SIGFPE";
    break;
  }

  // notify the user which signal was caught. We use printf, because this is the most
  // basic output function. Once you get a crash, it is possible that more complex output
  // systems like streams and the like may be corrupted. So we make the most basic call
  // possible to the lowest level, most standard print function.
  if (name)
    fprintf(stderr, "Caught signal %d (%s)\n", signum, name);
  else
    fprintf(stderr, "Caught signal %d\n", signum);

  Exception::PrintStackTrace();

  exit(signum);
}

void SignalHandler::SetupSignalHandlers()
{
  if (signal((int)SIGTERM, SignalHandler::ExitSignalHandler) == SIG_ERR)
    Exception::error(__FILE__, __LINE__,
                     std::string{"Error setting up signal handlers."},
                     static_cast<int>(ERROR_NUMBER::IGNORE));

  if (signal((int)SIGSEGV, SignalHandler::ExitSignalHandler) == SIG_ERR)
    Exception::error(__FILE__, __LINE__,
                     std::string{"Error setting up signal handlers."},
                     static_cast<int>(ERROR_NUMBER::IGNORE));

  if (signal((int)SIGINT, SignalHandler::ExitSignalHandler) == SIG_ERR)
    Exception::error(__FILE__, __LINE__,
                     std::string{"Error setting up signal handlers."},
                     static_cast<int>(ERROR_NUMBER::IGNORE));

  if (signal((int)SIGILL, SignalHandler::ExitSignalHandler) == SIG_ERR)
    Exception::error(__FILE__, __LINE__,
                     std::string{"Error setting up signal handlers."},
                     static_cast<int>(ERROR_NUMBER::IGNORE));

  if (signal((int)SIGABRT, SignalHandler::ExitSignalHandler) == SIG_ERR)
    Exception::error(__FILE__, __LINE__,
                     std::string{"Error setting up signal handlers."},
                     static_cast<int>(ERROR_NUMBER::IGNORE));

  if (signal((int)SIGFPE, SignalHandler::ExitSignalHandler) == SIG_ERR)
    Exception::error(__FILE__, __LINE__,
                     std::string{"Error setting up signal handlers."},
                     static_cast<int>(ERROR_NUMBER::IGNORE));
}
} // namespace DTDP