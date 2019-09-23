/* src/Communicator.cpp
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

#include "Communicator.hpp"

#include <string>

namespace DTDP
{
Communicator::Communicator(int argc, char **argv)
    : Communicator()
{
  aInitialized = static_cast<bool>(IsMPIinitialized());
  if (!aInitialized)
  {
    MPI_Init(&argc, &argv);
    aComm = MPI_COMM_WORLD;
    MPI_Comm_size(aComm, &aSize);
    MPI_Comm_rank(aComm, &aRank);
    aMasterRank = MASTER_RANK;
    aOwner = true;
  }
  else
    Exception::error(__FILE__, __LINE__,
                     std::string{"MPI execution environment initialized."});
}

Communicator::Communicator(const MPI_Comm &comm)
    : aComm{MPI_COMM_NULL},
      aSize{0},
      aRank{0},
      aMasterRank{MASTER_RANK},
      aInitialized{false},
      aOwner{false}
{
  if (IsMPIinitialized())
  {
    if (comm != MPI_COMM_NULL)
    {
      aInitialized = true;
      aComm = comm;
      MPI_Comm_size(aComm, &aSize);
      MPI_Comm_rank(aComm, &aRank);
    }
    // else // default to world comm
    // {
    //   aComm = MPI_COMM_WORLD;
    //   MPI_Comm_size(aComm, &aSize);
    //   MPI_Comm_rank(aComm, &aRank);
    //   aInitialized = true;
    // }
  }
}

Communicator::Communicator(const Communicator &c)
    : aComm{c.aComm},
      aSize{c.aSize},
      aRank{c.aRank},
      aMasterRank{c.aMasterRank},
      aInitialized{c.aInitialized},
      aOwner{false} {}

Communicator &Communicator::operator=(const Communicator &c)
{
  if (this != &c)
  {
    aComm = c.aComm;
    aSize = c.aSize;
    aRank = c.aRank;
    aMasterRank = c.aMasterRank;
    aInitialized = c.aInitialized;

    aOwner = false;
  }
  return *this;
}

Communicator::~Communicator()
{
  Finalize();
}

} // namespace DTDP
