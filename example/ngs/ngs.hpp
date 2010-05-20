// Copyright (C) 2008 - 2010 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

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
#include <cstring>
#include <algorithm>
#include <signal.h>

#ifndef ALPS_NGS_HPP
#define ALPS_NGS_HPP

namespace alps {
	class mcoptions {
		public:
			mcoptions(int argc, char* argv[]) : valid(false) {
				boost::program_options::options_description desc("Allowed options");
				desc.add_options()
					("help", "produce help message")
					("time-limit,T", boost::program_options::value<std::size_t>(&time_limit)->default_value(0), "time limit for the simulation")
					("max-bin-number,N", boost::program_options::value<std::size_t>(&max_bins)->default_value(0), "the maximum number of bins")
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
			std::size_t max_bins;
			std::size_t time_limit;
			std::string input_file;
	};
	// TODO: use boost::variant instead of string
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
				ar >> make_pvp("/parameters", *this);
			}
			mcparamvalue & operator[](std::string const & k) {
				return std::map<std::string, mcparamvalue>::operator[](k);
			}
			mcparamvalue const & operator[](std::string const & k) const {
				if (find(k) == end())
					throw std::invalid_argument("unknown argument: "  + k);
				return find(k)->second;
			}
			mcparamvalue value_or_default(std::string const & k, mcparamvalue const & v) const {
				if (find(k) == end())
					return mcparamvalue(v);
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
	template<typename T> class mcatomic {
		public:
			mcatomic() {}
			mcatomic(T const & v): value(v) {}
			mcatomic(mcatomic const & v): value(v.value) {}
			
			mcatomic & operator=(mcatomic const & v) {
				boost::lock_guard<boost::mutex> lock(mutex);
				value = v;
			}

			operator T() const { 
				boost::lock_guard<boost::mutex> lock(mutex);
				return value; 
			}
		private:
			T volatile value;
			boost::mutex mutable mutex;
	};
	template<typename Impl> class mcthreadsim : public Impl {
		public:
			mcthreadsim(typename Impl::parameters_type const & p)
				: Impl(p)
				, stop_flag(false)
			{}
			
			bool run(boost::function<bool ()> const & stop_callback) {
				boost::thread thread(boost::bind(&mcthreadsim<Impl>::thread_callback, this, stop_callback));
				Impl::run(boost::bind(&mcthreadsim<Impl>::run_callback, this));
				thread.join();
				return stop_flag;
			}
			
			virtual bool run_callback() {
				return stop_flag;
			}
			
			virtual void thread_callback(boost::function<bool ()> const & stop_callback) {
				while (true) {
					usleep(0.1 * 1e6);
					if (stop_flag = stop_callback())
						return;
				}
			}
		protected:
			mcatomic<bool> stop_flag;
			// boost::mutex mutex;
			// measurements and configuration need to be locked separately
	};
	template<typename Impl> class mcmpisim : public mcthreadsim<Impl> {
		public:
			enum {
				MPI_get_fraction	= 1,
				MPI_stop			= 2,
				MPI_collect			= 3,
				MPI_terminate		= 4
			};

			mcmpisim(typename mcthreadsim<Impl>::parameters_type const & p, boost::mpi::communicator const & c) 
				: mcthreadsim<Impl>(p)
				, next_check(8)
				, start_time(boost::posix_time::second_clock::local_time())
				, check_time(boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(next_check))
				, communicator(c)
			{
				MPI_Errhandler_set(communicator, MPI_ERRORS_RETURN);
			}

			double fraction_completed() {
				assert(communicator.rank() == 0);
				double fraction = mcthreadsim<Impl>::fraction_completed();
				int action = MPI_get_fraction;
				boost::mpi::broadcast(communicator, action, 0);
				boost::mpi::reduce(communicator, mcthreadsim<Impl>::fraction_completed(), fraction, std::plus<double>(), 0);
				return fraction;
			}

			typename mcthreadsim<Impl>::results_type collect_results() const {
				assert(communicator.rank() == 0);
				return collect_results(mcthreadsim<Impl>::result_names());
			}

			typename mcthreadsim<Impl>::results_type collect_results(typename mcthreadsim<Impl>::result_names_type const & names) const {
				int action = MPI_collect;
				char name[255];
				typename mcthreadsim<Impl>::results_type partial_results;
				// TODO ...
				return partial_results;
			}

			typename boost::shared_ptr<alea::mcdata<double> > collect_mcdata(std::string const & name) const {
				assert(name.size() < 255);
				int action = MPI_collect;
				boost::mpi::broadcast(communicator, action, 0);
				char buf[255];
				std::strcpy(buf, name.c_str());
				boost::mpi::broadcast(communicator, buf, 0);
				return boost::shared_ptr<alea::mcdata<double> >(reduce_results(name));
			}

			void terminate() const {
				if (communicator.rank() == 0) {
					int action = MPI_terminate;
					boost::mpi::broadcast(communicator, action, 0);
				}
			}

			void thread_callback(boost::function<bool ()> const & stop_callback) {
				if (communicator.rank() == 0) {
					for (bool flag = false; !flag && !stop_callback(); ) {
						usleep(0.1 * 1e6);
						if (!flag && boost::posix_time::second_clock::local_time() > check_time) {
							double fraction = fraction_completed();
							flag = (fraction >= 1.);
							next_check = std::min(2. * next_check, std::max(double(next_check), 0.8 * (boost::posix_time::second_clock::local_time() - start_time).total_seconds() / fraction * (1 - fraction)));
							check_time = boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(next_check);
							std::cerr << std::fixed << std::setprecision(1) << fraction * 100 << "% done next check in " << next_check << "s" << std::endl;
						}
					}
					mcthreadsim<Impl>::stop_flag = true;
					int action = MPI_stop;
					boost::mpi::broadcast(communicator, action, 0);
				} else
					process_requests();
			}

			void process_requests() {
				assert(communicator.rank() > 0);
				while (true) {
					int action;
					boost::mpi::broadcast(communicator, action, 0);
					switch (action) {
						case MPI_get_fraction:
							boost::mpi::reduce(communicator, mcthreadsim<Impl>::fraction_completed(), std::plus<double>(), 0);
							break;
						case MPI_collect:
							{
								char name[255];
								boost::mpi::broadcast(communicator, name, 0);
								reduce_results(name);
							}
							break;
						case MPI_stop:
						case MPI_terminate:
							mcthreadsim<Impl>::stop_flag = true;
							return;
					}
				}
			}
			
		private:
			alea::mcdata<double> * reduce_results(std::string const & name) const {
				using std::sqrt;
				using alps::numeric::sq;
				using alps::numeric::sqrt;
				alea::mcdata<double> data(dynamic_cast<AbstractSimpleObservable<double> const &>(
					mcthreadsim<Impl>::collect_results(typename result_names_type<mcthreadsim<Impl> >::type(1, name))[name]
				));
				data.set_bin_size(boost::mpi::all_reduce(communicator, data.bin_size(), boost::mpi::maximum<std::size_t>()));
				std::size_t min_bins = boost::mpi::all_reduce(communicator, data.bin_number(), boost::mpi::minimum<std::size_t>());
				uint64_t count, count_all;
				double mean, mean_all, error, error_all, variance, variance_all, binvalue;
				std::vector<double> binvalues;
				boost::tie(count, mean, error, variance, binvalue) = data.get_reduceable_data(boost::mpi::all_reduce(communicator, data.bin_number(), boost::mpi::minimum<std::size_t>()));
				if (communicator.rank()) {
					boost::mpi::reduce(communicator, count, std::plus<double>(), 0);
					boost::mpi::reduce(communicator, mean, std::plus<double>(), 0);
					boost::mpi::reduce(communicator, error, std::plus<double>(), 0);
					boost::mpi::reduce(communicator, variance, std::plus<double>(), 0);
					boost::mpi::gather(communicator, binvalue, 0);
					return NULL;
				} else {
					boost::mpi::reduce(communicator, count, count_all, std::plus<double>(), 0);
					boost::mpi::reduce(communicator, mean, mean_all, std::plus<double>(), 0);
					boost::mpi::reduce(communicator, error, error_all, std::plus<double>(), 0);
					boost::mpi::reduce(communicator, variance, variance_all, std::plus<double>(), 0);
					boost::mpi::gather(communicator, binvalue, binvalues, 0);
					return new alea::mcdata<double>(count_all, mean_all / double(count_all), sqrt(error_all) / double(count_all), variance_all / double(count_all), data.bin_size(), binvalues);
				}
			}
			
			int next_check;
			boost::posix_time::ptime start_time;
			boost::posix_time::ptime check_time;
			boost::mpi::communicator communicator;
	};
}

#endif