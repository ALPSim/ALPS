//  (C) Copyright 2010 Lukas Gamper <gamperl -at- gmail.com>
//  Use, modification, and distribution are subject to the Boost Software 
//  License, Version 1.0. (See at <http://www.boost.org/LICENSE_1_0.txt>.)

#ifndef ALPS_NGS_HPP
#define ALPS_NGS_HPP

#include <alps/alea.h>
#include <alps/hdf5.hpp>

#include <boost/bind.hpp>
#include <boost/timer.hpp>
#include <boost/utility.hpp>
#include <boost/function.hpp>
#include <boost/type_traits.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

#include <vector>
#include <unistd.h>
#include <signal.h>

#include <mpi.h>

namespace alps {
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
	template<typename S> double fraction_completed(S const & s) {
		return s.fraction_completed();
	}
	class mcoptions {
		public:
			mcoptions(int argc, char* argv[])
				: valid(false)
				, mpi(false)
			{
				boost::program_options::options_description desc("Allowed options");
				desc.add_options()
					("help", "produce help message")
					("mpi", "run in parallel using MPI") 
					("time-limit,T", boost::program_options::value<std::size_t>(&limit)->default_value(0),"time limit for the simulation")
					("input-file", boost::program_options::value<std::string>(&file), "input file");
				boost::program_options::positional_options_description p;
				p.add("input-file", 1);
				boost::program_options::variables_map vm;
				boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
				boost::program_options::notify(vm);
				valid = !vm.count("help");
				if (vm.count("help")) {
					std::cout << desc << std::endl;
					return;
				}
				if (vm.count("mpi"))
					mpi = true;
				if (file.empty())
					boost::throw_exception(std::runtime_error("No job file specified"));
			}
			bool is_valid() const { return valid; }
			int time_limit() const { return limit; }
			bool use_mpi() const { return mpi; }
			std::string input_file() const { return file; }
		private:
			bool valid;
			bool mpi;
			std::string file;
			std::size_t limit;
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
				return boost::lexical_cast<T>(*this);
			}
	};
	class mcparams : public std::map<std::string, mcparamvalue> {
		public: 
			mcparams(mcoptions const & o) {
				hdf5::iarchive ar(o.input_file());
				ar >> make_pvp("/parameters", this);
				operator[]("time_limit") = o.time_limit();
				operator[]("input_file") = o.input_file();
			}
			mcparamvalue & operator[](std::string const & k) {
				return std::map<std::string, mcparamvalue>::operator[](k);
			}
			mcparamvalue const & operator[](std::string const & k) const {
				if (find(k) == end())
					throw std::invalid_argument("unknown argument" + k);
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
	class mcbase {
		public:
			typedef mcparams parameters_type;
			typedef ObservableSet results_type;
			typedef std::vector<std::string> result_names_type;
			mcbase(parameters_type const & params): params(params) {}
			void save(boost::filesystem::path const & path) const {
				std::cout << "write checkpoint: " << path.file_string() << std::endl;
				boost::filesystem::path backup = boost::filesystem::exists(path) ? path.parent_path() / ( path.filename() + ".bak" ) : path;
				if (boost::filesystem::exists(backup))
					boost::filesystem::remove(backup);
				{
					hdf5::oarchive ar(backup.file_string());
					ar 
						<< make_pvp("/parameters", params)
						<< make_pvp("/simulation/realizations/0/clones/0/results", results);
				}
				if (backup != path) {
					boost::filesystem::remove(path);
					boost::filesystem::rename(backup, path);
				}
			}
			void load(boost::filesystem::path const & path) { 
				hdf5::iarchive ar(path.file_string());
				ar 
					>> make_pvp("/parameters", params)
					>> make_pvp("/simulation/realizations/0/clones/0/results", results);
			}
			result_names_type result_names() {
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
					partial_results[*it] = results[*it];
				return partial_results;
			}
		protected:
			parameters_type params;
			results_type results;
	};
	template <typename Impl> class mcrun {
		public:
			typedef typename Impl::parameters_type parameters_type;
			typedef typename Impl::results_type results_type;
			typedef typename Impl::result_names_type result_names_type;
			mcrun(parameters_type const & params, int argc = 0, char *argv[] = NULL)
				: impl(params)
				, finalized(false)
				, input_file(params["input_file"])
				, end_time(boost::posix_time::microsec_clock::local_time() + boost::posix_time::seconds(std::size_t(params["time_limit"])))
			{
				struct sigaction * action = new struct sigaction[3];
				for (std::size_t i = 0; i < 3; ++i) {
					memset(&action[i], 0, sizeof(action[i]));
					action[i].sa_handler = &mcrun::signal_hanler;
				}
				sigaction(SIGINT, &action[0], NULL);
				sigaction(SIGTERM, &action[1], NULL);
				sigaction(SIGXCPU, &action[2], NULL);
			}
			void save(boost::filesystem::path const & path) const { impl.save(path); }
			void load(boost::filesystem::path const & path) { impl.load(path); }
			bool run(boost::function<bool ()> const & stop_callback) {
				do {
					impl.do_update();
					impl.do_measurements();
				} while(!stop_callback() && fraction_completed() < 1);
				return fraction_completed() >= 1;
			}
			bool stop() {
				if ((signal || boost::posix_time::second_clock::local_time() > end_time) && !finalized) {
					save(checkpoint_name());
					finalized = true;
				}
				return finalized;
			}
			result_names_type result_names() const { return impl.result_names(); }
			result_names_type unsaved_result_names() const { return impl.unsaved_result_names(); }
			results_type collect_results() const { impl.collect_results(); }
			results_type collect_results(result_names_type const & names) const { return impl.collect_results(names); }
			double fraction_completed() const { return impl.fraction_completed(); }
		protected:
			virtual std::string checkpoint_name() const {
				return "sim.h5";
			}
		private:
			static void signal_hanler(int code) {
				std::cerr << "caught signal: " << code << std::endl;
				signal = true;
			}
			Impl impl;
			bool finalized;
			static bool signal;
			std::string input_file;
			boost::posix_time::ptime end_time;
	};
	template<typename Impl> bool mcrun<Impl>::signal = false;
	template <typename Impl> class mcmpirun : public mcrun<Impl> {
		public:
			mcmpirun(typename mcrun<Impl>::parameters_type const & params, int argc, char *argv[])
				: mcrun<Impl>(params)
			{
				MPI_Init (&argc, &argv);
				MPI_Comm_rank (MPI_COMM_WORLD, &rank);
				MPI_Comm_size (MPI_COMM_WORLD, &size);
				std::cerr << "starting simulation: " << rank << std::endl;
			}
			bool stop() {
				
				
				
				
				
//	s.save("sim.h5");
				
				
				
				
				
				return mcrun<Impl>::stop();
			}
		protected:
			std::string checkpoint_name() const {
				return "sim-" + boost::lexical_cast<std::string>(rank) + ".h5";
			}
		private:
			int rank;
			int size;
	};
}
#endif
