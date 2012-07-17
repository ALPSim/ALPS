In header ``#include <alps/ngs/hash.hpp>``.

The ALPS hash functions are an extention of  ``boost::hash``

Hash
----

.. doxygenstruct:: alps::hash


Example:

.. code-block:: c++

   alps::hash<double> hasher;
   std::size_t h = hasher(3.141); 

Hash combine
------------

.. doxygenfunction:: alps::hash_combine

Effect:

.. code-block:: c++

   s ^= hash_value(v) + 0x6dbc79f65d57ddf9ULL + (s << 12) + (s >> 4);

.. note::

   hash_value is call unqualified 
