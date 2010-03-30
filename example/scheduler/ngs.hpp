//  (C) Copyright 2010 Lukas Gamper <gamperl -at- gmail.com>
//  Use, modification, and distribution are subject to the Boost Software 
//  License, Version 1.0. (See at <http://www.boost.org/LICENSE_1_0.txt>.)

#ifndef ALPS_NGS_HPP
#define ALPS_NGS_HPP

#include <alps/alea.h>
#include <alps/hdf5.hpp>

#include <boost/utility.hpp>
#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/type_traits.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>
#include <algorithm>
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
			void reset(std::vector<char> const & buf) { 
				buffer = buf;
				pos = 0;
			}
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
			std::size_t size() const { return buffer.size(); }
			char & front() { return buffer.front(); }
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
					("tree-base,b", boost::program_options::value<std::size_t>(&base)->default_value(16), "number of children per node in the mpi communicationtree")
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
			std::size_t time_limit() const { return limit; }
			std::size_t tree_base() const { return base; }
			bool use_mpi() const { return mpi; }
			bool is_verbose() const { return verbose; }
			std::string input_file() const { return file; }
		private:
			bool valid;
			bool mpi;
			bool verbose;
			std::string file;
			std::size_t limit;
			std::size_t base;
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
				operator[]("tree_base") = o.tree_base();
				operator[]("input_file") = o.input_file();
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
#ifdef ALPS_HAVE_MPI
	enum mctags {
		MC_send_params		= 0x01,
		MC_get_fraction		= 0x02,
		MC_finalize			= 0x03
	};
	class mccommunicator {
		public:
			mccommunicator(
				  boost::function<void(int, boost::any &)> const & begin_
				, boost::function<void(int, boost::any &, mcidump &)> const & merge_
				, boost::function<void(int, boost::any &, mcodump &)> const & pack_
				, boost::function<void(int, boost::any &)> const & finish_
				, int argc
				, char *argv[]
				, int base = 16
			)
				: begin(begin_), merge(merge_), pack(pack_), finish(finish_)
			{
				MPI_Init (&argc, &argv);
				MPI_Comm_rank (MPI_COMM_WORLD, &rank);
				MPI_Comm_size (MPI_COMM_WORLD, &size);
				if(size > 1) {
					int layer = base;
					while (rank >= layer - 1)
						layer *= base;
					layer /= base;
					parent = std::max(0, (rank - layer + 1) / base + layer / base - 1);
					for (int i = 1; !rank && i < std::min(base - 1, size); ++i)
						children.push_back(i);
					for (int i = 0; i < base; ++i)
						if ((rank - layer + 1) * base + layer * base - 1 + i < size)
							children.push_back((rank - layer + 1) * base + layer * base - 1 + i);
				}
				std::stringstream out;
				out << "node: " << std::setw(4) << rank << ", parent: " << std::setw(4) << parent;
				if (children.size()) {
					out << ", children:";
					for (std::vector<int>::const_iterator it = children.begin(); it != children.end(); ++it)
						out << " " << std::setw(4) << *it;
				}
				out << std::endl;
				std::cerr << out.str();
			}
			bool is_master() {
				return !rank;
			}
			int get_rank() {
				return rank;
			}
			void broadcast(int tag, mcodump & buf) {
				if (rank)
					throw std::runtime_error("Only rank 0 can send a broadcast");
				tasks.resize(tasks.size() + 1);
				tasks.back().get<2>() = clock();
				begin(tag, tasks.back().get<1>());
				tag += (tasks.size() - 1) << 16;
				if (children.size()) {
					std::vector<MPI_Status> states(children.size()); 
					std::vector<MPI_Request> requests(children.size());
					int count = buf.size();
					for (std::size_t i = 0; i < children.size(); ++i)
						MPI_Isend(&(buf.front()), count, MPI_BYTE, children[i], tag, MPI_COMM_WORLD, &(requests[i]));
					MPI_Waitall(requests.size(), &(requests.front()), &(states.front()));
				}
			}
			void probe() {
				int flag;
				MPI_Status status;
				if (MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status) != MPI_SUCCESS)
					throw std::runtime_error("Error when receiving MPI message: " + boost::lexical_cast<std::string>(status.MPI_ERROR));
				if (flag) {
					int count;
					int tag = status.MPI_TAG & 0xFF;
					int index = status.MPI_TAG >> 16;
					MPI_Get_count(&status, MPI_BYTE, &count);
					std::vector<char> buf(count);
					MPI_Recv(&(buf.front()), count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
					if (status.MPI_SOURCE == parent) {
						tasks.resize(std::max(tasks.size(), std::size_t(index + 1)));
						begin(tag, tasks[index].get<1>());
						std::vector<MPI_Status> states(children.size()); 
						std::vector<MPI_Request> requests(children.size());
						for (std::size_t i = 0; i < children.size(); ++i)
							MPI_Isend(&(buf.front()), count, MPI_BYTE, children[i], status.MPI_TAG, MPI_COMM_WORLD, &(requests[i]));
						MPI_Waitall(requests.size(), &(requests.front()), &(states.front()));
					} else {
						if (std::find(children.begin(), children.end(), status.MPI_SOURCE) == children.end())
							throw std::runtime_error("Invalid MPI source: " + boost::lexical_cast<std::string>(status.MPI_SOURCE));
						mcidump ibuf;
						ibuf.reset(buf);
						merge(tag, tasks[index].get<1>(), ibuf);
						tasks[index].get<0>()++;
					}
					if (tasks[index].get<0>() == children.size()) {
						if (!rank) {
							std::cerr << "broadcast: " << std::setw(2) << tag << ", took: " << std::setprecision(4) << (clock() - tasks[index].get<2>()) / (double)CLOCKS_PER_SEC << "s" << std::endl;
							finish(tag, tasks[index].get<1>());
						} else if (tasks[index].get<0>() == children.size()) {
							mcodump obuf;
							pack(tag, tasks[index].get<1>(), obuf);
							count = obuf.size();
							MPI_Request request;
							MPI_Isend(&(obuf.front()), count, MPI_BYTE, parent, status.MPI_TAG, MPI_COMM_WORLD, &request);
							MPI_Wait(&request, &status);
						}
						tasks[index] = boost::make_tuple(0, boost::any(), 0);
					}
				}
			}
		private:
			int rank;
			int size;
			int parent;
			std::vector<int> children;
			std::vector<boost::tuple<int, boost::any, clock_t> > tasks;
			boost::function<void(int, boost::any &)> begin;
			boost::function<void(int, boost::any &)> finish;
			boost::function<void(int, boost::any &, mcidump &)> merge;
			boost::function<void(int, boost::any &, mcodump &)> pack;
	};
	template <typename Impl> class mcmpirun : public mcrun<Impl> {
		public:
			mcmpirun(typename mcrun<Impl>::parameters_type const & params, int argc, char *argv[])
				: mcrun<Impl>(params)
				, next_check(8)
				, checkpointing(false)
				, communicator(
					  boost::lambda::bind(&mcmpirun<Impl>::begin, boost::ref(*this), boost::lambda::_1, boost::lambda::_2)
					, boost::lambda::bind(&mcmpirun<Impl>::merge, boost::ref(*this), boost::lambda::_1, boost::lambda::_2, boost::lambda::_3)
					, boost::lambda::bind(&mcmpirun<Impl>::pack, boost::ref(*this), boost::lambda::_1, boost::lambda::_2, boost::lambda::_3)
					, boost::lambda::bind(&mcmpirun<Impl>::finish, boost::ref(*this), boost::lambda::_1, boost::lambda::_2)
					, argc
					, argv
					, params["tree_base"]
				)
				, start_time(boost::posix_time::second_clock::local_time())
				, check_time(boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(next_check))
			{}
			bool stop() {
				if (communicator.is_master() && next_check > 0 && boost::posix_time::second_clock::local_time() > check_time && !mcrun<Impl>::finalized) {
					next_check *= -1;
					mcodump dump;
					communicator.broadcast(MC_get_fraction, dump);
				}
				communicator.probe();
				return !checkpointing && mcrun<Impl>::stop();
			}
			void begin(int tag, boost::any & data) {
				switch (tag) {
					case MC_get_fraction:
						data = mcrun<Impl>::fraction_completed();
						break;
					case MC_finalize:
						data = std::vector<std::pair<std::size_t, typename mcrun<Impl>::results_type> >(1, make_pair(1, collect_results(*this)));
						mcrun<Impl>::save("sim-" + boost::lexical_cast<std::string>(communicator.get_rank()) + ".out.h5");
						break;
				}
			}
			void merge(int tag, boost::any & data, mcidump & buf) {
				switch (tag) {
					case MC_get_fraction:
						double fraction;
						buf >> fraction;
						data = fraction + boost::any_cast<double>(data);
						break;
					case MC_finalize:
						typename mcrun<Impl>::results_type result;
						buf >> result;
						std::vector<std::pair<std::size_t, typename mcrun<Impl>::results_type> > results = boost::any_cast<std::vector<std::pair<std::size_t, typename mcrun<Impl>::results_type> > >(data);
						results.push_back(make_pair(1, result));
						while (results.size() > 1 && results.back().first == (results.rbegin() + 1)->first) {
							(results.rbegin() + 1)->first *= 2;
							(results.rbegin() + 1)->second << results.back().second;
							results.pop_back();
						}
						data = results;
						break;
				}
			}
			void pack(int tag, boost::any & data, mcodump & buf) {
				switch (tag) {
					case MC_get_fraction:
						buf << boost::any_cast<double>(data);
						break;
					case MC_finalize:
						std::vector<std::pair<std::size_t, typename mcrun<Impl>::results_type> > results = boost::any_cast<std::vector<std::pair<std::size_t, typename mcrun<Impl>::results_type> > >(data);
						for (std::size_t i = results.size() - 1; i > 0; --i)
							results[i - 1].second << results[i].second;
						buf << results[0].second;
						checkpointing = false;
						break;
				}
			}
			void finish(int tag, boost::any & data) {
				switch (tag) {
					case MC_get_fraction:
						next_check = 8 + 0.5 * (boost::posix_time::second_clock::local_time() - start_time).total_seconds() / boost::any_cast<double>(data) * (1 - boost::any_cast<double>(data));
						std::cerr << std::fixed << std::setprecision(1) << boost::any_cast<double>(data) * 100 << "% done next check in " << next_check << "s" << std::endl;
						check_time = boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(next_check);
						if (boost::any_cast<double>(data) >= 1)
							finalize();
						break;
					case MC_finalize:
						std::vector<std::pair<std::size_t, typename mcrun<Impl>::results_type> > results = boost::any_cast<std::vector<std::pair<std::size_t, typename mcrun<Impl>::results_type> > >(data);
						for (std::size_t i = results.size() - 1; i > 0; --i)
							results[i - 1].second << results[i].second;
						boost::filesystem::path path = "sim.out.h5";
						if (mcrun<Impl>::verbose)
							std::cerr << "write file: " << path.file_string() << std::endl;
						boost::filesystem::path backup = boost::filesystem::exists(path) ? path.parent_path() / ( path.filename() + ".bak" ) : path;
						if (boost::filesystem::exists(backup))
							boost::filesystem::remove(backup);
						{
							hdf5::oarchive ar(backup.file_string());
							ar << make_pvp("/simulation/results", results[0].second);
						}
						if (backup != path) {
							boost::filesystem::remove(path);
							boost::filesystem::rename(backup, path);
						}
						checkpointing = false;
						break;
				}
			}
		protected:
			void finalize() {
				if (communicator.is_master()) {
					mcodump dump;
					communicator.broadcast(MC_finalize, dump);
				}
				mcrun<Impl>::finalized = true;
				checkpointing = true;
			}
		private:
			int next_check;
			bool checkpointing;
			mccommunicator communicator;
			boost::posix_time::ptime start_time;
			boost::posix_time::ptime check_time;
	};
#endif
}
#endif
