.. highlight:: c
ALPS Matrix class
=================

alps::numeric::matrix is a matrix class that models the concept of a [wiki:concepts/resizable_matrix Resizable Matrix] and can store various types T as its elements. For blas types (e.g. float) the general_matrix can be used directly for BLAS or LAPACK calls offered by boost::numeric::bindings.

In header ``#include <alps/numeric/matrix.hpp>``


Iterating over elements
-----------------------

If you want to iterate over the elements of the matrix, you have to be aware that the matrix is **column-major**,
which means that the the columns will be continuous in memory (having a unit stride).

.. important::
 Hence, the *innermost loop* should increment the *row* to allow the CPU to access the elements in an efficient way.

The preferred way of iterating over elements is using col_element_iterators, which iterate through the elements of a column (incrementing the row)::

 typedef std::pair<col_element_iterator,col_element_iterator> col_range_t;
 for(size_type j = 0; j < num_cols(m); ++j)
    for( col_range_t range = col(m,j) ; range.first != range.second ; ++range.first)
        do_something(*range.first);

Alternatively you may use element-wise access::

 for(size_type j = 0; j < num_cols(m); ++j)
    for( size_type i=0; i < num_rows(m) ++i)
        do_something(m(i,j));



Description
-----------

.. doxygenclass:: alps::numeric::matrix
   :members: 



