
#include <alps/alea.h>
#include <alps/hdf5.hpp>

#include <boost/mpi.hpp>
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/utility.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>

#include <map>
#include <vector>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <signal.h>

namespace alps {
	class mcidump : public IDump {
		public:
			mcidump(std::vector<char> const & buf) : buffer(buf), pos(0) {}
			#define ALPS_MCIDUMP_DO_TYPE(T)																		\
				void read_simple(T & x) { read_buffer(&x, sizeof(T)); }											\
				void read_array(std::size_t n, T * p) { read_buffer(p, n * sizeof(T)); }
			ALPS_MCIDUMP_DO_TYPE(bool)
			ALPS_MCIDUMP_DO_TYPE(char)
			ALPS_MCIDUMP_DO_TYPE(signed char)
			ALPS_MCIDUMP_DO_TYPE(unsigned char)
			ALPS_MCIDUMP_DO_TYPE(short)
			ALPS_MCIDUMP_DO_TYPE(unsigned short)
			ALPS_MCIDUMP_DO_TYPE(int)
			ALPS_MCIDUMP_DO_TYPE(unsigned int)
			ALPS_MCIDUMP_DO_TYPE(long)
			ALPS_MCIDUMP_DO_TYPE(unsigned long)
			#ifdef BOOST_HAS_LONG_LONG
				ALPS_MCIDUMP_DO_TYPE(long long)
				ALPS_MCIDUMP_DO_TYPE(unsigned long long)
			#endif
			ALPS_MCIDUMP_DO_TYPE(float)
			ALPS_MCIDUMP_DO_TYPE(double)
			ALPS_MCIDUMP_DO_TYPE(long double)
			# undef ALPS_MCIDUMP_DO_TYPE
			void read_string(std::size_t n, char * x) { read_buffer(x, n); }
		private:
			void read_buffer(void * p, std::size_t n) {
				if(pos + n > buffer.size())
					throw std::runtime_error("read past buffer");
				std::memcpy(p, &(buffer[pos]), n);
				pos += n;
			}
			std::size_t pos;
			std::vector<char> buffer;
	};
	class mcodump : public ODump {
		public:
			std::vector<char> const & data() { return buffer; }
			#define ALPS_MCODUMP_DO_TYPE(T)																		\
				void write_simple(T x) { write_buffer(&x, sizeof(T)); }											\
				void write_array(std::size_t n, T const * p) { write_buffer(p, n * sizeof(T)); }
			ALPS_MCODUMP_DO_TYPE(bool)
			ALPS_MCODUMP_DO_TYPE(char)
			ALPS_MCODUMP_DO_TYPE(signed char)
			ALPS_MCODUMP_DO_TYPE(unsigned char)
			ALPS_MCODUMP_DO_TYPE(short)
			ALPS_MCODUMP_DO_TYPE(unsigned short)
			ALPS_MCODUMP_DO_TYPE(int)
			ALPS_MCODUMP_DO_TYPE(unsigned int)
			ALPS_MCODUMP_DO_TYPE(long)
			ALPS_MCODUMP_DO_TYPE(unsigned long)
			#ifdef BOOST_HAS_LONG_LONG
				ALPS_MCODUMP_DO_TYPE(long long)
				ALPS_MCODUMP_DO_TYPE(unsigned long long)
			#endif
			ALPS_MCODUMP_DO_TYPE(float)
			ALPS_MCODUMP_DO_TYPE(double)
			ALPS_MCODUMP_DO_TYPE(long double)
			# undef ALPS_MCODUMP_DO_TYPE
			void write_string(std::size_t n, const char * x) { write_buffer(x, n); }
		private:
			void write_buffer(void const * p, std::size_t n) {
				std::size_t write_pos = buffer.size();
				buffer.resize(write_pos + n);
				std::memcpy(&(buffer[write_pos]), p, n);
			}
			std::vector<char> buffer;
	};
	class mcoptions {
		public:
			mcoptions(int argc, char* argv[]) : valid(false) {
				boost::program_options::options_description desc("Allowed options");
				desc.add_options()
					("help", "produce help message")
					("time-limit,T", boost::program_options::value<std::size_t>(&time_limit)->default_value(0), "time limit for the simulation")
					("input-file", boost::program_options::value<std::string>(&input_file), "input file in hdf5 format");
				boost::program_options::positional_options_description p;
				p.add("input-file", 1);
				boost::program_options::variables_map vm;
				boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
				boost::program_options::notify(vm);
				if (!(valid = !vm.count("help")))
					std::cout << desc << std::endl;
				else if (input_file.empty())
					boost::throw_exception(std::runtime_error("No job file specified"));
			}
			bool valid;
			std::string input_file;
			std::size_t time_limit;
	};
	class mcparamvalue : public std::string {
		public:
			mcparamvalue() {}
			template <typename T> mcparamvalue(mcparamvalue const & v)
				: std::string(v) 
			{}
			template <typename T> mcparamvalue(T const & v)
				: std::string(boost::lexical_cast<std::string>(v)) 
			{}
			template <typename T> typename boost::disable_if<typename boost::is_base_of<std::string, T>::type, mcparamvalue &>::type operator=(T const & v) {
				std::string::operator=(boost::lexical_cast<std::string>(v));
				return *this;
			}
			template <typename T> typename boost::enable_if<typename boost::is_base_of<std::string, T>::type, mcparamvalue &>::type operator=(T const & v) {
				std::string::operator=(v);
				return *this;
			}
			template <typename T> operator T() const {
				return boost::lexical_cast<T, std::string>(*this);
			}
	};
	class mcparams : public std::map<std::string, mcparamvalue> {
		public: 
			mcparams(std::string const & input_file) {
				hdf5::iarchive ar(input_file);
				ar >> make_pvp("/parameters", this);
			}
			mcparamvalue & operator[](std::string const & k) {
				return std::map<std::string, mcparamvalue>::operator[](k);
			}
			mcparamvalue const & operator[](std::string const & k) const {
				if (find(k) == end())
					throw std::invalid_argument("unknown argument: "  + k);
				return find(k)->second;
			}
			void serialize(hdf5::oarchive & ar) const {
				for (const_iterator it = begin(); it != end(); ++it)
					ar << make_pvp(it->first, static_cast<std::string>(it->second));
			}
			void serialize(hdf5::iarchive & ar) {
				std::vector<std::string> list = ar.list_children(ar.get_context());
				for (std::vector<std::string>::const_iterator it = list.begin(); it != list.end(); ++it) {
					std::string v;
					ar >> make_pvp(*it, v);
					insert(std::make_pair(*it, v));
				}
			}
	};
	class mcsignal{
		public:
			mcsignal() {
				static bool initialized;
				if (!initialized) {
					static struct sigaction action;
					initialized = true;
					memset(&action, 0, sizeof(action));
					action.sa_handler = &mcsignal::slot;
					sigaction(SIGINT, &action, NULL);
					sigaction(SIGTERM, &action, NULL);
					sigaction(SIGXCPU, &action, NULL);
					sigaction(SIGQUIT, &action, NULL);
					sigaction(SIGUSR1, &action, NULL);
					sigaction(SIGUSR2, &action, NULL);
				}
			}
			operator bool() const { return which; }
			static int code() { return which; }
			static void slot(int signal) { which = signal; }
		private:
			static int which;
	};
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	template<typename S> struct result_names_type {
		typedef typename S::result_names_type type;
	};
	template<typename S> struct results_type {
		typedef typename S::results_type type;
	};
	template<typename S> struct parameters_type {
		typedef typename S::parameters_type type;
	};
	template<typename S> typename result_names_type<S>::type result_names(S const & s) {
		return s.result_names();
	}
	template<typename S> typename result_names_type<S>::type unsaved_result_names(S const & s) {
		return s.unsaved_result_names();
	}
	template<typename S> typename results_type<S>::type collect_results(S const & s) {
		return s.collect_results();
	}
	template<typename S> typename results_type<S>::type collect_results(S const & s, typename result_names_type<S>::type const & names) {
		return s.collect_results(names);
	}
	template<typename S> typename results_type<S>::type collect_results(S const & s, std::string const & name) {
		return collect_results(s, typename result_names_type<S>::type(1, name));
	}
	template<typename S> double fraction_completed(S const & s) {
		return s.fraction_completed();
	}
	class mcbase {
		public:
			typedef mcparams parameters_type;
			typedef ObservableSet results_type;
			typedef std::vector<std::string> result_names_type;
			mcbase(parameters_type const & p)
				: params(p)
				, random(boost::mt19937(), boost::uniform_real<>())
			{}
			virtual void do_update() = 0;
			virtual void do_measurements() = 0;
			virtual double fraction_completed() const = 0;
			void save(boost::filesystem::path const & path) const {
				boost::filesystem::path original = path.parent_path() / (path.filename() + ".h5");
				boost::filesystem::path backup = path.parent_path() / (path.filename() + ".bak");
				if (boost::filesystem::exists(backup))
					boost::filesystem::remove(backup);
				{
					hdf5::oarchive ar(backup.file_string());
					ar 
						<< make_pvp("/parameters", params)
						<< make_pvp("/simulation/realizations/0/clones/0/results", results);
				}
				if (boost::filesystem::exists(original))
					boost::filesystem::remove(original);
				boost::filesystem::rename(backup, original);
			}
			void load(boost::filesystem::path const & path) { 
				hdf5::iarchive ar(path.file_string());
				ar >> make_pvp("/simulation/realizations/0/clones/0/results", results);
			}
			bool run(boost::function<bool ()> const & stop_callback) {
				do {
					do_update();
					do_measurements();
				} while(!stop_callback() && fraction_completed() < 1);
				return fraction_completed() >= 1;
			}
			result_names_type result_names() const {
				result_names_type names;
				for(results_type::const_iterator it = results.begin(); it != results.end(); ++it)
					names.push_back(it->first);
				return names;
			}
			result_names_type unsaved_result_names() const { return result_names_type(); }
			results_type collect_results() const { return results; }
			results_type collect_results(result_names_type const & names) const {
				results_type partial_results;
				for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it)
					partial_results << results[*it];
				return partial_results;
			}
		protected:
			parameters_type params;
			results_type results;
			boost::variate_generator<boost::mt19937, boost::uniform_real<> > random;
	};

// write atomic class
// template <class T> class atomic { atomic(), operator=, operator T ; private: T volatile, mutex}

	template<typename Impl> class mcthreadsim : public Impl {
		public:
			mcthreadsim(typename mcwrapsim<Impl>::parameters_type const & p)
				: mcwrapsim<Impl>(p)
				, stop_flag(false)
			{}
			bool run(boost::function<bool ()> const & stop_callback) {
				boost::thread thread(boost::bind(&mcthreadsim<Impl>::thread_callback, this, stop_callback));
				mcwrapsim<Impl>::run(boost::bind(&mcthreadsim<Impl>::run_callback, this));
				thread.join();
				return stop_flag;
			}
			virtual bool run_callback() {
				boost::lock_guard<boost::mutex> lock(mutex);
				return stop_flag;
			}
			virtual void thread_callback(boost::function<bool ()> const & stop_callback) {
				while (true) {
					usleep(0.1 * 1e6);
					{
                        bool new_stop_flag = stop_callback();
                        {
                            boost::lock_guard<boost::mutex> lock(mutex);
                            stop_flag = new_stop_flag;
                        }
						if (new_stop_flag)
							return;
					}
				}
			}
		protected:
			bool volatile stop_flag;
			boost::mutex mutex;
            // measurements and configuration need to be locked separately
	};
	template<typename Result> class mcmpimerge {
		public:
			std::vector<char> operator()(std::vector<char> const & arg1, std::vector<char> const & arg2) {
				mcidump idump1(arg1), idump2(arg2);
				Result results1, results2;
				idump1 >> results1;
				idump2 >> results2;
				results1 << results2;
				mcodump odump;
				odump << results1;
				return odump.data();
			}
	};
}
namespace boost { 
	namespace mpi {
		template<typename Result> struct is_commutative<alps::mcmpimerge<Result>, std::vector<char> > : public mpl::true_ { };
	} 
}
namespace alps {
	template<typename Impl> class mcmpisim : public mcthreadsim<Impl> {
		public:
			enum {
				MPI_get_fraction		= 1,
				MPI_stop				= 2,
				MPI_collect				= 3,
				MPI_terminate			= 4
			};
            
			mcmpisim(typename mcthreadsim<Impl>::parameters_type const & p, boost::mpi::communicator const & c) 
				: mcthreadsim<Impl>(p)
				, next_check(8)
				, start_time(boost::posix_time::second_clock::local_time())
				, check_time(boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(next_check))
				, communicator(c)
			{}
            
			double fraction_completed() {
				assert (communicator.rank()==0);
                
                double fraction = mcthreadsim<Impl>::fraction_completed();
                int action = MPI_get_fraction;
                boost::mpi::broadcast(communicator, action, 0);
                boost::mpi::reduce(communicator, mcthreadsim<Impl>::fraction_completed(), fraction, std::plus<double>(), 0);
                return fraction;
			}
            
			typename mcthreadsim<Impl>::results_type collect_results() const {
				assert (communicator.rank()==0);

                int action = MPI_collect;
                std::vector<char> buf;
                mcodump odump;
                odump << mcthreadsim<Impl>::collect_results();
                boost::mpi::broadcast(communicator, action, 0);
                boost::mpi::reduce(communicator, odump.data(), buf, mcmpimerge<typename mcthreadsim<Impl>::results_type>(), 0);
                mcidump idump(buf);
                typename mcthreadsim<Impl>::results_type results;
                idump >> results;
                return results;
			}
            
			typename mcthreadsim<Impl>::results_type collect_results(typename mcthreadsim<Impl>::result_names_type const & names) const { 
				typename mcthreadsim<Impl>::results_type results = collect_results(), partial_results;
				for(typename mcthreadsim<Impl>::result_names_type::const_iterator it = names.begin(); it != names.end(); ++it)
					partial_results << results[*it];
				return partial_results;
			}
            
			void terminate() const {
				if (communicator.rank()==0) {
					int action = MPI_terminate;
					boost::mpi::broadcast(communicator, action, 0);
				}
			}
            
			void thread_callback(boost::function<bool ()> const & stop_callback) {
				if (communicator.rank()==0) {
					bool flag = false;
					while (!flag && !stop_callback()) {
						usleep(0.1 * 1e6);
						if (!flag && boost::posix_time::second_clock::local_time() > check_time) {
							double fraction = fraction_completed();
							flag = (fraction >= 1.);
							next_check = std::min(2. * next_check, std::max(double(next_check), 0.8 * (boost::posix_time::second_clock::local_time() - start_time).total_seconds() / fraction * (1 - fraction)));
							check_time = boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(next_check);
							std::cerr << std::fixed << std::setprecision(1) << fraction * 100 << "% done next check in " << next_check << "s" << std::endl;
						}
					}
					{
						boost::lock_guard<boost::mutex> lock(mcthreadsim<Impl>::mutex);
						flag = mcthreadsim<Impl>::stop_flag = true;
					}
					int action = MPI_stop;
					boost::mpi::broadcast(communicator, action, 0);
				} else
					process_requests();
			}
            
			void process_requests() {
				assert (communicator.rank()>0);
                
				while (true) {
					int action;
					boost::mpi::broadcast(communicator, action, 0);
					switch (action) {
						case MPI_get_fraction:
							boost::mpi::reduce(communicator, mcthreadsim<Impl>::fraction_completed(), std::plus<double>(), 0);
							break;
						case MPI_collect:
							{
								mcodump odump;
								odump << collect_results();
								boost::mpi::reduce(communicator, odump.data(), mcmpimerge<typename mcthreadsim<Impl>::results_type>(), 0);
							}
							break;
						case MPI_stop:
						case MPI_terminate:
							{
								boost::lock_guard<boost::mutex> lock(mcthreadsim<Impl>::mutex);
								mcthreadsim<Impl>::stop_flag = false;
								return;
							}
					}
				}
			}
		private:
			int next_check;
			boost::posix_time::ptime start_time;
			boost::posix_time::ptime check_time;
			boost::mpi::communicator communicator;
	};
}
