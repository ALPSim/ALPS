// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <string>
#include <vector>
#include <sstream>
#ifndef TRACE
#define TRACE
namespace mocasito {
	#define MOCASITO_TRACE																	\
		static char const * __trace_file__ = __FILE__;										\
		static std::string const __trace_function__(__FUNCTION__);							\
		mocasito::trace<> __trace__ = mocasito::trace<>(__LINE__, __trace_file__, __trace_function__.c_str());
	#define MOCASITO_IO_THROW(message)														\
		{																					\
			MOCASITO_TRACE;																	\
			throw(std::runtime_error(message + mocasito::trace<>::str()));					\
		}
	template<typename T = void> class trace {
		public:
			trace(std::size_t line, char const * file, char const * function)
				: _line(line), _file(file), _function(function)
			{
				_stack.push_back(this);
			}
			~trace() {
				_stack.pop_back();
			}
			static std::string str() {
				std::size_t cnt = 0;
				std::ostringstream buffer;
				buffer << std::endl;
				for (typename std::vector<trace<T> *>::reverse_iterator it = _stack.rbegin(); it != _stack.rend(); ++it, ++cnt) {
					std::string file = (**it)._file;
					buffer << "#" << cnt << " " << file.substr(file.find_last_of('/') == std::string::npos ? 0 : file.find_last_of('/')) << " line " << (**it)._line << " in " << (**it)._function << std::endl;
				}
				return buffer.str();
			}
		private:
			std::size_t const _line;
			char const * _file;
			char const * _function;
			static std::vector<trace<T> *> _stack;
	};
	template<typename T> std::vector<trace<T> *> trace<T>::_stack;
	
}
#endif
