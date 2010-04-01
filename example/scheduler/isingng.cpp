# include "isingng.hpp"

//#define SINGLE_RUN
//#define THREAD_RUN
#define MPI_RUN

#ifdef SINGLE_RUN
	typedef ising_sim sim_type;
#endif
#ifdef THREAD_RUN
	typedef alps::mcthreadsim<ising_sim> sim_type;
#endif
#ifdef MPI_RUN
	typedef alps::mcmpisim<ising_sim> sim_type;
#endif
bool stop_callback(boost::posix_time::ptime const & end_time) {
	return alps::mcsignal() || boost::posix_time::second_clock::local_time() > end_time;
}
int main(int argc, char *argv[]) {
	alps::mcoptions options(argc, argv);
	if (options.valid) {
		sim_type::parameters_type params(options.input_file);
#ifdef MPI_RUN
		boost::mpi::environment env(argc, argv);
		boost::mpi::communicator c;
		sim_type s(params, c);
#else
		sim_type s(params);
#endif
		s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(options.time_limit)));
#ifdef MPI_RUN
		s.save("sim-" + boost::lexical_cast<std::string>(c.rank()));
		if (s.is_master()) {
			sim_type::results_type results = collect_results(s);
			{
				alps::hdf5::oarchive ar("sim.h5");
				ar << make_pvp("/parameters", params) << make_pvp("/simulation/results", results);
			}
			s.terminate();
		} else
			s.process_requests();
#else
		s.save("sim");
		std::cout << collect_results(s, "Magnetization");
#endif
	}
}
