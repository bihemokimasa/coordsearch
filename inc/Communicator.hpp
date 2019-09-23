/* inc/Communicator.hpp
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

// provides basic routines for an MPI (Message Passing Interface) execution environment

#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <mpi.h>

#include "DynProg.hpp"
#include "Exception.hpp"

namespace DTDP
{
class Communicator
{
public:
  Communicator(const Communicator &c);
  Communicator &operator=(const Communicator &c);
  ~Communicator();

  Communicator(int argc, char **argv);
  Communicator(const MPI_Comm &comm = MPI_COMM_WORLD);
  Index Size() const { return aSize; }
  Index Rank() const { return aRank; }
  static Index Size(const MPI_Comm &comm)
  {
    int sz;
    MPI_Comm_size(comm, &sz);
    return sz;
  }

  static Index Rank(const MPI_Comm &comm)
  {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
  }

  Index MasterRank() const { return aMasterRank; }

  bool Initialized() const { return aInitialized; }
  operator MPI_Comm() const { return aComm; }

  static bool IsMPIinitialized()
  {
    Index initialized;
    MPI_Initialized(&initialized);
    return initialized;
  }

  static bool IsMPIfinalized()
  {
    Index flag;
    MPI_Finalized(&flag);
    return flag;
  }

  Index Abort() const
  {
    Index err{0};
    MPI_Abort(aComm, err);
    return err;
  }

  void Finalize()
  {
    if (!IsMPIfinalized() && aOwner)
      MPI_Finalize();
  }

  bool SeqOrMaster() const // sequential program or master process
  {
    return !aInitialized || (aInitialized && aRank == aMasterRank);
  }

private:
  MPI_Comm aComm; // communication context
  Index aSize;
  Index aRank;
  Index aMasterRank;
  bool aInitialized;
  bool aOwner;
};
} // namespace DTDP

#endif
