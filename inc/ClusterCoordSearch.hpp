/* inc/ClusterCoordSearch.hpp
*
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


#ifndef CLUSTERCOORDSEARCH_HPP
#define CLUSTERCOORDSEARCH_HPP

#include <vector>

#include "Communicator.hpp"
#include "GPSCoordSearch.hpp"

namespace DTDP
{
class ClusterCoordSearch : virtual public GPSCoordSearch
{
public:
  ClusterCoordSearch() = default;
  ClusterCoordSearch(const ClusterCoordSearch &) = default;
  ClusterCoordSearch &operator=(const ClusterCoordSearch &) = default;
  virtual ~ClusterCoordSearch();

  explicit ClusterCoordSearch(const GPSCoordSearch &gcs,
                              Index num_teams,
                              const MPI_Comm &wc);

  void SetCommunicator(Index num_teams, const MPI_Comm &wc);

  virtual void init_fval() override;
  virtual void poll() override;

  virtual void lhs_search() override;
  virtual void update_pattern() override { return; /* DON'T pattern.sort();*/ };
  virtual void print_pattern(std::ostream &os = std::cout) const override;
  virtual void print_progress(std::ostream &os) override;
  virtual void print_result(std::ostream &os) override;

  bool seek_best_nearby(Index &fcalls_local, Index &fcalls);

  void reset_recv_buffer();

  /* column count if ptfxnew_col_count */
  Index ptfxnew_col_count() const { return {get_n() + 1}; }

private:
  /* workspace receive buffers */
  std::vector<D> ptfxnew; /* workspace variable:  we will store candidate trial points */

  Index team_count; /* number of teams */
  Index team_id;    /* rank of team */
  // std::vector<Index> coords; /* in the cartesian topology for the MPI environment * /

  Communicator world_comm;
  MPI_Comm mpi_world_comm; /* world communicator */

  Communicator sub_comm;
  MPI_Comm mpi_sub_comm; /* specific to a team */
};
} // namespace DTDP

#endif
