// SPDX-License-Identifier: MIT
#include <alps/utility/cli.hpp>
#include <alps/utility/citations.hpp>
#include <alps/utility/copyright.hpp>
#include <iostream>
#include <stdexcept>
#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif

alps::cli_mpi_guard::cli_mpi_guard(int& argc, char**& argv) : owns_mpi_(false) {
#ifdef ALPS_HAVE_MPI
  int initialized = 0, finalized = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (finalized)
    throw std::runtime_error("Cannot run a CLI query after MPI has finalized");
  if (!initialized) {
    MPI_Init(&argc, &argv);
    owns_mpi_ = true;
  }
#else
  (void)argc;
  (void)argv;
#endif
}

alps::cli_mpi_guard::~cli_mpi_guard() {
#ifdef ALPS_HAVE_MPI
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (owns_mpi_ && !finalized) MPI_Finalize();
#else
  (void)owns_mpi_;
#endif
}

bool alps::cli_is_master() {
#ifdef ALPS_HAVE_MPI
  int initialized = 0, finalized = 0, rank = 0;
  MPI_Initialized(&initialized);
  MPI_Finalized(&finalized);
  if (initialized && !finalized) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
#else
  return true;
#endif
}

bool alps::handle_cli_information(int argc, char** argv, const std::string& component,
                                const std::function<void()>& print_help) {
  std::string query;
  int queries = 0;
  bool extra_argument = false;
  for (int i = 1; i < argc; ++i) {
    const std::string argument(argv[i]);
    if (argument == "--") {
      extra_argument = true;
      break;
    }
    if (argument == "--help" || argument == "-h" || argument == "--license" ||
        argument == "-l" || argument == "--citations") {
      query = argument;
      ++queries;
    } else if (argument != "--mpi") {
      extra_argument = true;
    }
  }
  if (!queries) return false;
  if (queries != 1 || extra_argument)
    throw std::invalid_argument("Use one of --help, --license, or --citations without calculation arguments (optional --mpi)");

  cli_mpi_guard mpi(argc, argv);
  if (cli_is_master()) {
    if (query == "--help" || query == "-h") {
      print_help();
      std::cout << "\nInformation options (no simulation input required):\n"
                << "  --help, -h    show help\n"
                << "  --license, -l print license conditions\n"
                << "  --citations   print recommended citations\n";
    } else if (query == "--license" || query == "-l") {
      print_license(std::cout);
    } else {
      print_citation_details(std::cout, component);
    }
  }
  return true;
}
