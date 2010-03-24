//  (C) Copyright 2010 Lukas Gamper <gamperl -at- gmail.com>
//  Use, modification, and distribution are subject to the Boost Software 
//  License, Version 1.0. (See at <http://www.boost.org/LICENSE_1_0.txt>.)

/* $Id: abstract_task.C 3822 2010-01-30 22:02:39Z troyer $ */

#ifndef ALPS_NGS_HPP
#define ALPS_NGS_HPP

#include <alps/alea.h>
#include <alps/hdf5.hpp>
#include <alps/parameter.h>

#include <boost/function.hpp>
#include <boost/filesystem/path.hpp>

#include <vector>

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
	template <typename Impl> class mcrun {
		public:
			typedef Parameters parameters_type;
			typedef ObservableSet results_type;
			typedef std::vector<std::string> result_names_type;
			mcrun(parameters_type const & parms): impl(parms) {}
			void save(boost::filesystem::path const & path) const { impl.save(path); }
			void load(boost::filesystem::path const & path) { impl.load(path); }
			bool run(boost::function<bool ()> const & stop_callback) {
				do {
					impl.do_update();
					impl.do_measurements();
				} while(!stop_callback() && fraction_completed() < 1);
				return fraction_completed() >= 1;
			}
			result_names_type result_names() const { return impl.result_names(); }
			result_names_type unsaved_result_names() const { return impl.unsaved_result_names(); }
			results_type collect_results() const { impl.collect_results(); }
			results_type collect_results(result_names_type const & names) const { return impl.collect_results(names); }
			double fraction_completed() const { return impl.fraction_completed(); }
		private:
			Impl impl;
	};
	class mcbase {
		public:
			typedef Parameters parameters_type;
			typedef ObservableSet results_type;
			typedef std::vector<std::string> result_names_type;
			mcbase(parameters_type const & parms): parms(parms) {}
			void save(boost::filesystem::path const & path) const {
				boost::filesystem::path backup = boost::filesystem::exists(path) ? path.parent_path() / ( path.filename() + ".bak" ) : path;
				if (boost::filesystem::exists(backup))
					boost::filesystem::remove(backup);
				{
					hdf5::oarchive ar(backup.file_string());
					ar 
						<< make_pvp("/parameters", parms)
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
					>> make_pvp("/parameters", parms)
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
			parameters_type parms;
			results_type results;
	};
}
#endif
