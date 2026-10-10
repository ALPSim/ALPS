rem
rem Copyright (C) 2009 - 2011 Matthias Troyer
rem
rem Distributed under the Boost Software License, Version 1.0. (See accompany-
rem ing file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

mkdir deplibs
mkdir deplibs\bin
mkdir deplibs\lib
mkdir deplibs\include
mkdir deplibs\share
mkdir deplibs\etc
copy ..\..\opt\lib deplibs\lib
copy ..\..\opt\bin deplibs\bin
copy ..\..\opt\include deplibs\include
