Signal Handler
==============
In header ``#include <alps/ngs/signal.hpp>``

The ``alps::ngs::signal`` class provides a nice interface to posix signals. On Windows no signals were captured.

``alps::ngs::signal`` is a singletone class, i.e. the signal handler are initialized the first time ``alps::ngs::signal`` is created. 



Description
-----------

.. doxygenclass:: alps::ngs::signal
   :members: 
   


