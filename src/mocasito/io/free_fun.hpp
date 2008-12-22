// Copyright (C) 2008 Lukas Gamper <gamperl -at- gmail.com>
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include "context.hpp"
#include "traits.hpp"
#include "util.hpp"
#ifndef IO_FREE_FUN
#define IO_FREE_FUN
namespace mocasito {
	namespace io {
		template<typename E, typename T> T & assign(T & v, context<E> const & c) {
			if (!c.is_attribute())
				set_size<T>::apply(v, c.extends()[0]);
			c.get(get_data<T>::apply(v));
			return v;
		}
		template<typename E, typename T> context<E> & assign(context<E> & c, T const & v) {
			if (has_fixed_size<T>::value && fixed_size<T>::value == 1)
				c.set(v);
			else
				c.set(get_data<T>::apply(v), get_size<T>::apply(v));
			return c;
		}
		template<typename E> context<E> & assign(context<E> & c, char const * v) {
			c.set(v, std::strlen(v));
			return c;
		}
		template<typename E, typename T> context<E> & append(context<E> & c, T const & v) { 
			c.append(get_data<T>::apply(v), get_size<T>::apply(v));
			return c;
		}
	}
}
#endif
