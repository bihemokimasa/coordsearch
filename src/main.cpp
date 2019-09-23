
#include <iostream>

#include "CoordSearch.hpp"
#include "GPSCoordSearch.hpp"
#include "testfunc.hpp"
#include "HJDirSearch.hpp"
#include "Communicator.hpp"
#include "ClusterCoordSearch.hpp"
#include "SignalHandler.hpp"

/* Paraboloid centered on (p[0],p[1]), with
   scale factors (p[2],p[3]) and minimum p[4] */

double my_f(const DTDP::Array1D<DTDP::D> &x,
            const DTDP::SolverParam &params)
{
  return Rosebrock(x.begin(), x.size());
}

int main(int argc, char *argv[])
{
  DTDP::Communicator comm(argc, argv);
  MPI_Comm wc = MPI_COMM_WORLD;

  DTDP::SignalHandler signalHandler;
  // Register signal handler to handle kill signal
  signalHandler.SetupSignalHandlers();

  // jb_main();
  const int n{2};
  double x[n] = {-1.2, 1.0};

  // DTDP::Index pattern_type{static_cast<DTDP::Index>(DTDP::CoordSearch::Pattern::Type::GPS2N)};
  DTDP::Index pattern_type{static_cast<DTDP::Index>(DTDP::CoordSearch::Pattern::Type::GPSNp1)};
  // DTDP::Index search_method{static_cast<DTDP::Index>(DTDP::GPSCoordSearch::SearchMethod::LINE_SEARCH)};
  DTDP::Index search_method{static_cast<DTDP::Index>(DTDP::GPSCoordSearch::SearchMethod::LHS_SEARCH)};
  // DTDP::Index search_method{static_cast<DTDP::Index>(DTDP::GPSCoordSearch::SearchMethod::NONE)};

  DTDP::Index lhs_pt_count{n + 1};

  // gcs.set_max_fcalls(8409);
  // gcs.set_max_iters(4981);

  if (DTDP::Communicator::Rank(wc) == 0)
  {

    DTDP::CoordSearch cs(&x[0], n);
    cs.set_init_slen(0.3);
    cs.search();

    std::cout << "Coordinate Search COMPLETE!" << std::endl
              << std::endl;

    DTDP::GPSCoordSearch gcs(&x[0],
                             n,
                             pattern_type,
                             1.0,
                             lhs_pt_count,
                             search_method);

    gcs.set_init_slen(0.3);
    gcs.search();

    std::cout << std::endl
              << "GPS Coordinate Search COMPLETE!" << std::endl
              << std::endl;

    DTDP::Array1D<double> startpt(&x[0], n), endpt(n);
    DTDP::HJDirSearch hj(startpt, &my_f);
    hj.SetRho(0.5);
    hj.SetDisplay(true);
    hj.SetMaxFunEvals(200);
    hj.hooke(DTDP::SolverParam());

    std::cout << "HOOKE and JEEVES Search COMPLETE!" << std::endl
              << std::endl;
  }

  try
  {
    DTDP::GPSCoordSearch gcs(&x[0],
                             n,
                             pattern_type,
                             1.0,
                             lhs_pt_count,
                             search_method);

    gcs.set_max_fcalls(200);
    gcs.set_init_slen(0.3);
    DTDP::ClusterCoordSearch parcs(gcs, 2, wc);

    parcs.search();
  }
  catch (const DTDP::Exception &e)
  {
    if (comm.Initialized())
    {
      std::cerr << "Rank: "
                << comm.Rank()
                << " "
                << e.what()
                << std::endl;
      DTDP::Exception::PrintStackTrace();
      comm.Abort();
    }
    else
    {
      std::cerr << e.what() << std::endl;
      DTDP::Exception::PrintStackTrace();
      exit(1);
    }
  }
  // parcs.test();

  return 0;
}
