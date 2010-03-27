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

#include <iomanip>
#include <vector>
//#include <unistd.h>
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
	class mcidump : public IDump {
		public:
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
			std::size_t size() { return buffer.size(); }
			char * data() { return &buffer.front(); }
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
			std::size_t size() { return buffer.size(); }
			char * data() { return &buffer.front(); }
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
			mcoptions(int argc, char* argv[])
				: valid(false)
				, mpi(false)
				, verbose(false)
			{
				boost::program_options::options_description desc("Allowed options");
				desc.add_options()
					("help", "produce help message")
					("mpi", "run in parallel using MPI") 
					("time-limit,T", boost::program_options::value<std::size_t>(&limit)->default_value(0), "time limit for the simulation")
					("verbose,v", "verbose mode")
					("input-file", boost::program_options::value<std::string>(&file), "input file in hdf5 format");
				boost::program_options::positional_options_description p;
				p.add("input-file", 1);
				boost::program_options::variables_map vm;
				boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
				boost::program_options::notify(vm);
				valid = !vm.count("help");
				if (vm.count("mpi"))
					mpi = true;
				if (vm.count("verbose"))
					verbose = true;
				if (vm.count("help"))
					std::cout << desc << std::endl;
				else if (file.empty())
					boost::throw_exception(std::runtime_error("No job file specified"));
			}
			bool is_valid() const { return valid; }
			int time_limit() const { return limit; }
			bool use_mpi() const { return mpi; }
			bool is_verbose() const { return verbose; }
			std::string input_file() const { return file; }
		private:
			bool valid;
			bool mpi;
			bool verbose;
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
				return boost::lexical_cast<T, std::string>(*this);
			}
	};
	class mcparams : public std::map<std::string, mcparamvalue> {
		public: 
			mcparams(mcoptions const & o) {
				hdf5::iarchive ar(o.input_file());
				ar >> make_pvp("/parameters", this);
				operator[]("time_limit") = o.time_limit();
				operator[]("verbose") = o.is_verbose();
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
				if (bool(params["verbose"]))
					std::cerr << "write file: " << path.file_string() << std::endl;
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
				, verbose(params["verbose"])
				, input_file(params["input_file"])
				, end_time(boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(std::size_t(params["time_limit"])))
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
			void save(boost::filesystem::path const & path) const { 
				impl.save(path); 
			}
			void load(boost::filesystem::path const & path) { impl.load(path); }
			bool run(boost::function<bool ()> const & stop_callback) {
				do {
					impl.do_update();
					impl.do_measurements();
				} while(!stop_callback() && fraction_completed() < 1);
				return fraction_completed() >= 1;
			}
			bool stop() {
				if ((signal || (boost::posix_time::second_clock::local_time() > end_time)) && !finalized)
					finalize();
				return finalized;
			}
			result_names_type result_names() const { return impl.result_names(); }
			result_names_type unsaved_result_names() const { return impl.unsaved_result_names(); }
			results_type collect_results() const { return impl.collect_results(); }
			results_type collect_results(result_names_type const & names) const { return impl.collect_results(names); }
			double fraction_completed() const { return impl.fraction_completed(); }
		protected:
			virtual void finalize() {
				if (mcrun<Impl>::verbose)
					std::cerr << "checkpoint after: " << std::fixed << std::setprecision(1) << fraction_completed() * 100 << "%" << std::endl;
				save(input_file + ".out.h5");
				finalized = true;
			}
			bool finalized;
			bool verbose;
		private:
			static void signal_hanler(int code) {
				std::cerr << "caught signal: " << code << std::endl;
				signal = true;
			}
			Impl impl;
			static bool signal;
			std::string input_file;
			boost::posix_time::ptime end_time;
	};
	template<typename Impl> bool mcrun<Impl>::signal = false;
	enum mctags {
		MC_send_params		= 0x01,
		MC_get_fraction		= 0x02,
		MC_finalize			= 0x03
	};
	template <typename Impl> class mcmpirun : public mcrun<Impl> {
		public:
			mcmpirun(typename mcrun<Impl>::parameters_type const & params, int argc, char *argv[])
				: mcrun<Impl>(params)
				, source(0)
				, next_check(3)
				, check_time(boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(next_check))
			{
				MPI_Init (&argc, &argv);
				MPI_Comm_rank (MPI_COMM_WORLD, &rank);
				MPI_Comm_size (MPI_COMM_WORLD, &size);
				int nmesg = 2; // size > 16 ? 16 : size < 4 ? 2 : size / 2;
				if(size > 1) {
					int nstep = nmesg;
					if (rank > 0)
						while((rank % nstep) == 0)
							nstep *= nstep;
						source = (rank / nstep) * nstep;
					nstep = 1;
					while (nstep < size)
						nstep *= nmesg;
					do {
						if (rank % nstep == 0)
							for (int i = 1; i < nmesg; ++i)
								if (rank + i * (nstep / nmesg) < size)
									targets.push_back(rank + i * (nstep / nmesg));
						nstep /= nmesg;
					} while (nstep != 1);
				}
				if (mcrun<Impl>::verbose) {
					std::cerr << "start on rank: " << std::setw(3) << rank << ", source: " << std::setw(3) << source;
					if (targets.size()) {
						std::cerr << ", targets:";
						for (std::vector<int>::const_iterator it = targets.begin(); it != targets.end(); ++it)
							std::cerr << " " << std::setw(3) << *it;
					}
					std::cerr << std::endl;
				}
			}
			bool stop() {
				if (!rank && boost::posix_time::second_clock::local_time() > check_time && !mcrun<Impl>::finalized) {
					check_time = boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(next_check);
					
					
					
					if (mcrun<Impl>::verbose)
						std::cerr << "checkin progress" << std::endl;
			
				}
				return mcrun<Impl>::stop();
			}
		protected:
			void finalize() {
				if (!rank) {
					mcodump dump;
					dump << collect_results(*this);
					
//					std::cout << dump.size() << std::endl;
					
				}
				mcrun<Impl>::save("sim-" + boost::lexical_cast<std::string>(rank) + ".out.h5");
				mcrun<Impl>::finalized = true;
			}
		private:
/*			void broadcast(int tag, mcodump const & obuf, std::vector<mcidump> & ibufs) {
				MPI_Status status;
				int flag = true;
				if (rank != 0){
					MPI_Iprobe(source, tag, MPI_COMM_WORLD, &flag, &status);
					if (flag) {
						int count;
						MPI_Get_count(&status, MPI_BYTE, &count);
						ibuf.resize(count);
						MPI_Recv(&(buf[0]), count, MPI_BYTE, source, status.MPI_TAG, MPI_COMM_WORLD, &status);
						for (int k = 0; k < targets.size(); ++k) {
							MPI_Request * request = new MPI_Request; 
							MPI_Isend(&(buf[0]), count, MPI_BYTE, targets[k], status.MPI_TAG, MPI_COMM_WORLD, request);
							
							
							MPI_Request_free(request);
							
							
						}
					}
				} else {
					for(int k = 0; k < targets.size(); ++k) {
						MPI_Request * request = new MPI_Request; 
						int count = buf.size();
						MPI_Isend(&(buf[0]), count, MPI_BYTE, targets[k], tag, MPI_COMM_WORLD, request);
						
						
						MPI_Request_free(request);
						
						
					}
				}
				return "";
			}
*/			int rank;
			int size;
			int source;
			int next_check;
			std::vector<int> targets;
			boost::posix_time::ptime check_time;
			
	};
}
#endif
