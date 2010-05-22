# include "ising.hpp"

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
void dump_result(std::string const & name, alps::alea::mcdata<double> const & result, std::string const & rank = "") {
    std::stringstream s;
    s << std::fixed << std::setprecision(5) 
      << name << rank << ": "
      << " count: " << result.count()
      << " mean: " << result.mean()
      << " error: " << result.error()
      << " variance: " << result.variance()
      << " tau: " << result.tau()
      << " bins: [" << result.bins().front() << ",...," << result.bins().back() << "]"
      << " #bins: " << result.bin_number()
      << " bin size: " << result.bin_size()
      << std::endl;
    std::cout << s.str();
}
void dump_result(std::string const & name, alps::alea::mcdata<std::vector<double> > const & result, std::string const & rank = "") {
    std::stringstream s;
    s << std::fixed << std::setprecision (5) 
      << name << rank << ": "
      << " count: " << result.count()
      << " mean: [" << result.mean().front() << ",...," << result.mean().back() << "]"
      << " error: [" << result.error().front() << ",...," << result.error().back() << "]"
      << " variance: [" << result.variance().front() << ",...," << result.variance().back() << "]"
      << " tau: [" << result.tau().front() << ",...," << result.tau().back() << "]"
      << " #bins: " << result.bin_number()
      << " bin size: " << result.bin_size()
      << std::endl;
    std::cout << s.str();
}
int main(int argc, char *argv[]) {
    alps::mcoptions options(argc, argv);
    if (options.valid) {
        sim_type::parameters_type params(options.input_file);
        params["maxbinnumber"] = options.max_bins;
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
        alps::results_type<sim_type>::type results = s.collect_local_results();
        for (alps::results_type<sim_type>::type::const_iterator it = results.begin(); it != results.end(); ++it)
            if (dynamic_cast<alps::alea::mcdata<double> const *>(it->second) != NULL)
                dump_result(it->first, *dynamic_cast<alps::alea::mcdata<double> const *>(it->second), "(" + boost::lexical_cast<std::string>(c.rank()) + ")");
            else if (dynamic_cast<alps::alea::mcdata<std::vector<double> > const *>(it->second) != NULL)
                dump_result(it->first, *dynamic_cast<alps::alea::mcdata<std::vector<double> > const *>(it->second), "rank: " + boost::lexical_cast<std::string>(c.rank()));
            else
                boost::throw_exception(std::runtime_error("unknown observable type"));
        if (c.rank()==0) {
            {
                alps::results_type<sim_type>::type results = collect_results(s);
                alps::hdf5::oarchive ar("sim.h5");
                ar << alps::make_pvp("/parameters", params);
                for (alps::results_type<sim_type>::type::const_iterator it = results.begin(); it != results.end(); ++it)
                    ar << alps::make_pvp("/simulation/results/" + it->first, *(it->second));
                for (alps::results_type<sim_type>::type::const_iterator it = results.begin(); it != results.end(); ++it)
                    if (dynamic_cast<alps::alea::mcdata<double> const *>(it->second) != NULL)
                        dump_result(it->first, *dynamic_cast<alps::alea::mcdata<double> const *>(it->second));
                    else if (dynamic_cast<alps::alea::mcdata<std::vector<double> > const *>(it->second) != NULL)
                        dump_result(it->first, *dynamic_cast<alps::alea::mcdata<std::vector<double> > const *>(it->second));
                    else
                        boost::throw_exception(std::runtime_error("unknown observable type"));
            }
            s.terminate();
        } else
            s.process_requests();
#else
        s.save("sim");
        collect_results(s);
        std::cout << collect_results(s, "Magnetization");
#endif
    }
}
