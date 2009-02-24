// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <string>
#include <vector>
#include <sstream>
#ifndef IO_UTIL
#define IO_UTIL
namespace mocasito {
	namespace io {
		#define MOCASITO_TRACE mocasito::io::trace<> __trace__ = mocasito::io::trace<>(__LINE__, __FUNCTION__, __FILE__);
		#define MOCASITO_IO_THROW(message)													\
			{																				\
				MOCASITO_TRACE;																\
				throw(std::runtime_error(message + mocasito::io::trace<>::str()));			\
			}
		namespace detail {
			struct trace_entry {
				trace_entry (std::size_t ln, std::string fn, std::string fl)
					: line(ln), function(fn), file(fl)
				{}
				std::size_t line;
				std::string function, file;
			};
		}
		template<typename T = void> class trace {
			public:
				trace(std::size_t line, std::string function, std::string file) {
					_stack.push_back(detail::trace_entry(line, function, file));
				}
				~trace() {
					_stack.pop_back();
				}
				static void push(std::size_t line, std::string function, std::string file) {
					_stack.push_back(detail::trace_entry(line, function, file));
				}
				static std::string str() {
					std::size_t cnt = 0;
					std::ostringstream buffer;
					for (std::vector<detail::trace_entry>::reverse_iterator it = _stack.rbegin(); it != _stack.rend(); ++it, ++cnt)
						buffer << "#" << cnt << " " << it->file << " line " << it->line << " in " << it->function << std::endl;
					return buffer.str();
				}
			private:
				static std::vector<detail::trace_entry> _stack;
		};
		template<typename T> std::vector<detail::trace_entry> trace<T>::_stack;
		#define MOCASITO_IO_FOREACH_SCALAR(callback)										\
			callback(char)																	\
			callback(signed char)															\
			callback(unsigned char)															\
			callback(short)																	\
			callback(unsigned short)														\
			callback(long)																	\
			callback(int)																	\
			callback(unsigned)																\
			callback(unsigned long)															\
			callback(long long)																\
			callback(unsigned long long)													\
			callback(float)																	\
			callback(double)																\
			callback(long double)															\
			callback(bool)
		namespace detail {
			typedef enum node_t { IO_UNKNOWN = -1, IO_GROUP, IO_DATA, IO_GROUP_ATTR, IO_DATA_ATTR } node_t;
		}
	}
}
#endif
