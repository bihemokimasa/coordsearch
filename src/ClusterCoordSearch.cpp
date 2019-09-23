/* src/ClusterCoordSearch.cpp
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

#include "ClusterCoordSearch.hpp"

namespace DTDP
{

ClusterCoordSearch::
    ClusterCoordSearch(const GPSCoordSearch &gcs,
                       Index num_teams,
                       const MPI_Comm &wc)
    : GPSCoordSearch{GPSCoordSearch(gcs)},
      ptfxnew{std::vector<D>(ptfxnew_col_count(), get_infinity())},
      team_count{num_teams},
      team_id{std::numeric_limits<int>::max()},
      world_comm{Communicator()},
      mpi_world_comm{MPI_COMM_NULL},
      sub_comm{Communicator()},
      mpi_sub_comm{MPI_COMM_NULL}
{
  const Index num_proc{Communicator::Size(wc)};
  const bool mpi_initialized{Communicator().IsMPIinitialized()};

  if (!mpi_initialized || num_proc == 0 || team_count < 1)
    Exception::error(__FILE__, __LINE__,
                     std::string{"MPI environment unitialized or np = 0 team_count < 2."});

  if (num_proc % (num_proc / num_teams) != 0)
    Exception::error(__FILE__, __LINE__,
                     std::string{"required equal np per team."});

  SetCommunicator(num_teams, wc);

  const Index maxdims{2}; /* required*/
  Index coords[maxdims];
  MPI_Cart_coords(mpi_world_comm, world_comm.Rank(), maxdims, &coords[0]);

  /* set team ID */
  team_id = coords[0];

  if (world_comm.Rank() == world_comm.MasterRank())
  {
    ptfxnew = std::vector<D>(team_count * ptfxnew_col_count(), get_infinity());
  }
}

void ClusterCoordSearch::reset_recv_buffer()
{
  /*
  *
  * reset the values ptfxnew from previous iteration
  * 
  */
  for (size_t i{0}; i < ptfxnew.size(); ++i)
  {
    ptfxnew[i] = get_infinity();
  }
}

ClusterCoordSearch::~ClusterCoordSearch()
{
  MPI_Comm_free(&mpi_sub_comm);
  MPI_Comm_free(&mpi_world_comm);
}

void ClusterCoordSearch::
    SetCommunicator(Index num_teams,
                    const MPI_Comm &wc)
{
  /*
  *
  * set teams of process by use of the MPI cartesian topology
  * where rows represent team ID
  * and 
  * columns represent the process belonging to the same team
  * 
  * we'll reorder the ranks 
  * and 
  * wrap around so that the grid is periodic
  * 
  */

  const Index num_proc{Communicator::Size(wc)};
  const Index cart_ndims{2}, rows{num_teams};
  const Index cols{num_proc / rows};

  const Index dim_sizes[cart_ndims] = {rows, cols};
  const Index wrap_around[cart_ndims] = {1, 1};
  const Index reorder{1};

  MPI_Cart_create(wc, cart_ndims, &dim_sizes[0],
                  &wrap_around[0], reorder,
                  &mpi_world_comm);

  world_comm = Communicator(mpi_world_comm);

  Index free_coords[cart_ndims] = {0, 1};

  MPI_Cart_sub(mpi_world_comm, &free_coords[0], &mpi_sub_comm);

  sub_comm = Communicator(mpi_sub_comm);
}

void ClusterCoordSearch::init_fval()
{
  /* initialize fx */
  if (get_search_method() == static_cast<Index>(SearchMethod::LHS_SEARCH))
  {
    /* take the lowest point from the poll around the initial pt 
    *
    * excludes the initial pt
    * 
    */ 
    lhs_search();
  }
  else
  {
    CoordSearch::init_fval();
  }
}

void ClusterCoordSearch::lhs_search()
{
  /* search based on the latin hypercube sampling 
  *
  * see poll() for documentation
  * 
  */

  reset_recv_buffer();

  Index i_row{team_id};

  Index fcalls_local{0}, fcalls{0};
  while (1)
  {
    GPSCoordSearch::lhs_search(i_row, fcalls_local);

    if (seek_best_nearby(
            fcalls_local,
            fcalls))
    {
      /* we received an abort signal since
      no function evaluation was made */
      break;
    }

    inc_fcalls(fcalls);
    accept_search_pt();

    // if (get_fxnew() < get_fx())
    // {
    //   set_fx(get_fxnew());
    //   set_pt(get_ptnew());
    // }

    i_row += team_count;
  }
}

void ClusterCoordSearch::poll()
{
  /* poll nearby points*/

  reset_recv_buffer();

  set_fxptbefore();

  /*
  * index to unit coordinate direction on the search pattern
  * 
  * a team operates on the direction coordinate whose ID is the team's ID
  */
  Index i_row{team_id};
  Index fcalls_local{0}, fcalls{0}; /* function evaluation counter */
  while (1)
  {
    // Index fcalls{0};
    GPSCoordSearch::poll(i_row, fcalls_local);

    if (seek_best_nearby(
            fcalls_local,
            fcalls))
    {
      /* we received an abort signal since
      no function evaluation was made */
      break;
    }

    inc_fcalls(fcalls);

    if (accept_ptnew())
    {
      break;
    }

    // if (get_fxnew() < get_fxbefore())
    // {
    //   update_fxptbefore();
    //   break;
    // }

    i_row += team_count;
  }

  update();
}

bool ClusterCoordSearch::
    seek_best_nearby(Index &fcalls_local,
                     Index &fcalls)
{
  /* we are going to share our new point values 
  and
  take the best new point

  to accumulate the total number of function calls 
  use MPI_Allreduce
  and set fcalls_local = 0 where process is non master 
  for the team.
  */

  if (team_count == 1)
  {
    fcalls = fcalls_local;
    if (fcalls == 0)
    {
      return true;
    }
    return false;
  }

  if (sub_comm.Rank() == sub_comm.MasterRank())
  {
    copy(get_ptnew(), ptfxnew);
    ptfxnew[get_n()] = get_fxnew();
  }
  else
  {
    fcalls_local = 0;
  }

  /* compute total of funcation calls */
  MPI_Allreduce(&fcalls_local, &fcalls,
                1, MPI_INT, MPI_SUM, mpi_world_comm);

  if (fcalls == 0)
  {
    return true;
  }

  const Index recv_count{ptfxnew_col_count()}; /* we goint to receive n + 1 points*/
  {
    /* global master process receives from other masters the newest points */
    Index coords[2] = {0, 0};                         /* required variable */
    if (world_comm.Rank() == world_comm.MasterRank()) /* coordinate (0,0) */
    {
      const Index req_count{team_count - 1}; /* number of requests to make */
      MPI_Request req[req_count];            /* required variable */
      MPI_Status stat[req_count];            /* required variable */

      for (Index t_id{1}; t_id < team_count; ++t_id)
      {
        Index src{0}, tag{t_id}; /* required variable */
        coords[0] = t_id;
        MPI_Cart_rank(mpi_world_comm, &coords[0], &src);

        MPI_Irecv(&ptfxnew.data()[t_id * recv_count],
                  recv_count, MPI_DOUBLE,
                  src, tag, mpi_world_comm, &req[t_id - 1]);
      }

      /* exit only after requests were completed */
      MPI_Waitall(req_count, &req[0], &stat[0]);
    }
    else if (sub_comm.Rank() == sub_comm.MasterRank()) /* coordinate (i, 0), i = 1, 2, ..., team_count -1 */
    {
      MPI_Request req;
      MPI_Status stat;

      MPI_Cart_coords(mpi_world_comm, world_comm.Rank(), 2, &coords[0]);
      Index dest{world_comm.MasterRank()}, tag{coords[0]};

      MPI_Isend(&ptfxnew.data()[0], ptfxnew.size(), MPI_DOUBLE,
                dest, tag, mpi_world_comm, &req);

      /* exit only after request was completed */
      MPI_Wait(&req, &stat);
    }
  }

  /* master process will not compute the minimum of the new points */
  if (world_comm.Rank() == world_comm.MasterRank())
  {
    for (Index t_id{1}; t_id < team_count; ++t_id)
    {
      D fxnewmin = ptfxnew[get_n() + t_id * recv_count]; /* temporary */

      if (fxnewmin < get_fxnew())
      {
        ptfxnew[get_n()] = ptfxnew[get_n() + t_id * recv_count];
        for (Index j{0}; j < get_n(); ++j)
          ptfxnew[j] = ptfxnew[j + t_id * recv_count];
      }
    }
  }

  MPI_Bcast(ptfxnew.data(), recv_count, MPI_DOUBLE,
            world_comm.MasterRank(), mpi_world_comm);

  set_fxnew(ptfxnew[get_n()]);
  set_ptnew(ptfxnew.data(), get_n());

  return false;
}

void ClusterCoordSearch::print_pattern(std::ostream &os) const
{
  if (world_comm.Rank() == world_comm.MasterRank())
  {
    GPSCoordSearch::print_pattern(os);
  }
}

void ClusterCoordSearch::print_progress(std::ostream &os)
{
  if (world_comm.Rank() == world_comm.MasterRank())
  {
    GPSCoordSearch::print_progress(os);
  }
}

void ClusterCoordSearch::print_result(std::ostream &os)
{
  if (world_comm.Rank() == world_comm.MasterRank())
  {
    GPSCoordSearch::print_result(os);
  }
}

// template <class T>
// void ClusterCoordSearch::
//     receive_from_other_masters(
//         const T *send_data,
//         Index send_count,
//         T *recv_data,
//         Index recv_count,
//         const MPI_Datatype &TYPE) const
// {

//   Index coords[2] = {0, 0};

//   if (world_comm.Rank() == world_comm.MasterRank()) /* coordinate (0,0) */
//   {
//     MPI_Request req[team_count - 1];
//     for (Index r{1}; r < team_count; r++)
//     {
//       Index src{0};
//       Index tag{r};
//       coords[0] = r;
//       MPI_Cart_rank(mpi_world_comm, &coords[0], &src);
//       MPI_Irecv(&recv_data[r * recv_count],
//                 recv_count,
//                 TYPE,
//                 src,
//                 tag,
//                 mpi_world_comm,
//                 &req[r - 1]);
//     }
//     MPI_Status stat[team_count - 1];
//     MPI_Waitall(team_count - 1, &req[0], &stat[0]);
//   }
//   else if (sub_comm.Rank() == sub_comm.MasterRank()) /* coordinate (i, 0), i = 1, 2, ..., team_count -1 */
//   {
//     MPI_Request req;

//     MPI_Cart_coords(mpi_world_comm, world_comm.Rank(), 2, &coords[0]);
//     Index dest{world_comm.MasterRank()};
//     Index tag{coords[0]};

//     MPI_Isend(&send_data[0],
//               send_count,
//               TYPE,
//               dest,
//               tag,
//               mpi_world_comm,
//               &req);

//     MPI_Status stat;
//     MPI_Wait(&req, &stat);
//   }
// }

// template <class T>
// void ClusterCoordSearch::
//     share_to_other_masters(const T *send_data,
//                            T *recv_data,
//                            Index count,
//                            const MPI_Datatype &TYPE) const
// {
//   if (world_comm.Rank() == world_comm.MasterRank())
//   {
//     MPI_Request req[team_count - 1];
//     Index coords[2] = {0, sub_comm.MasterRank()};
//     for (Index r{1}; r < team_count; ++r)
//     {
//       // determine rank for team r at coord(r, column==subcomm master rank)
//       coords[0] = r;
//       Index dest{0};
//       MPI_Cart_rank(mpi_world_comm, &coords[0], &dest);

//       MPI_Isend(&send_data[0],
//                 count,
//                 TYPE,
//                 dest,
//                 world_comm.Rank(),
//                 mpi_world_comm,
//                 &req[r - 1]);
//     }

//     MPI_Status stat[team_count - 1];
//     MPI_Waitall(team_count - 1, &req[0], &stat[0]);
//   }
//   else
//   {
//     if (sub_comm.Rank() == sub_comm.MasterRank()) /* coordinate (i, 0), i = 1, 2, ..., team_count-1 */
//     {

//       MPI_Request req;

//       MPI_Irecv(&recv_data[0],
//                 count,
//                 TYPE,
//                 world_comm.MasterRank(),
//                 world_comm.MasterRank(),
//                 mpi_world_comm,
//                 &req);

//       MPI_Status stat;
//       MPI_Wait(&req, &stat);
//     }
//   }
// }

// void ClusterCoordSearch::receive_ptfx_from_other_masters()
// {
//   // deduce number of function evaluations

//   Index reqs_count{0};

//   Index i_min{0};
//   Index coords[2] = {0, 0};

//   if (world_comm.Rank() == world_comm.MasterRank()) /* coordinate (0,0) */
//   {
//     gfx.data()[0] = get_fx();

//     MPI_Request fx_req[team_count - 1];
//     for (Index r{1}; r < team_count; r++)
//     {
//       if (moffset[r] < pattern_size())
//       {
//         Index src{0};
//         Index tag{r};
//         coords[0] = r;
//         MPI_Cart_rank(mpi_world_comm, &coords[0], &src);
//         MPI_Irecv(&gfx.data()[r],
//                   1,
//                   MPI_DOUBLE,
//                   src,
//                   tag,
//                   mpi_world_comm,
//                   &fx_req[r - 1]);

//         reqs_count++;
//       }
//     }
//     MPI_Status stat[reqs_count];
//     MPI_Waitall(reqs_count, &fx_req[0], &stat[0]);
//   }
//   else if (sub_comm.Rank() == sub_comm.MasterRank()) /* coordinate (i, 0), i = 1, 2, ..., aNumTeams-1 */
//   {
//     if (offset < pattern_size())
//     {
//       MPI_Request req[2];

//       MPI_Cart_coords(mpi_world_comm, world_comm.Rank(), 2, &coords[0]);
//       Index dest{world_comm.MasterRank()};
//       Index tag{coords[0]};
//       D fx{get_fx()};

//       MPI_Isend(&fx,
//                 1,
//                 MPI_DOUBLE,
//                 dest,
//                 tag,
//                 mpi_world_comm,
//                 &req[0]);
//     }
//   }

//   MPI_Bcast(&reqs_count,
//             1, MPI_INT,
//             world_comm.MasterRank(),
//             mpi_world_comm);

//   inc_fcalls(reqs_count);

//   MPI_Bcast(gfx.data(),
//             reqs_count,
//             MPI_DOUBLE,
//             world_comm.MasterRank(),
//             mpi_world_comm);
// }

} // namespace DTDP