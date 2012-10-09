#include <boost/concept_archetypes.hpp>

namespace alps {
namespace numeric {
namespace concepts {


template <typename T>
class matrix_archetype
{
  /**
     \brief Class docs?
   **/
private:
    typedef matrix_archetype self;
public:
    typedef T               value_type;     //< The type of the coefficients stored in the matrix
    typedef std::size_t     size_type;      //< An unsigned integral type used to represent the size (num_rows, num_cols) of the matrix
    typedef std::ptrdiff_t  difference_type;//< A signed integral type used to represent the distance between two of the container's iterators

    // TODO more typedefs

    /**
      * \brief Element access
      * @param i row index of the element
      * @param j column index of the element
      * Returns a reference to the element located at (i,j)
      *
      * @precond i < num_rows(m) && j < num_cols(m)
      * @postcond i < num_rows(m) && j < num_cols(m)
      *
      * @new_in{2.1}
      * We have introduced concept archetypes with automatic 
      * documentation derivation as well the ability to list new and 
      * changed things.
      **/
    value_type& operator()(size_type i, size_type j) { return value_type(); }

    /**
      * \brief Element access (constant)
      * @param i row index of the element
      * @param j column index of the element
      * @precond i < num_rows(m) && j < num_cols(m)
      * @new_in{2.3}
      * Returns a const reference to the element located at (i,j)
      **/
    value_type const& operator()(size_type i, size_type j) const { return value_type(); }

    /**
      * \brief Assignement operator
      * Assigns the matrix to the argument
      * @return A reference to this.
      *
      * @postcond The matrix has the same dimensions as m and the same coefficients.
      * @invariant m remains unchanged.
      *
      **/
    matrix_archetype& operator = (matrix_archetype const& m) { return *this; }

    /**
      * \brief Plus-assignemnt
      * Adds the matrix m to the matrix
      *
      * @precond The matrices have the same dimensions
      * @postcond TODO
      *
      * @return A reference to this.
      */
    matrix_archetype& operator += (matrix_archetype const& m){ return *this; }

    // TODO complete this
};

/**
  * \brief Returns the number of rows
  * @invariant m remains unchanged
 **/
template <typename T>
typename matrix_archetype<T>::size_type num_rows(matrix_archetype<T> const& m) { return typename matrix_archetype<T>::size_type(0); }

/**
  * \brief Returns the number of columns
  * @invariant m remains unchanged
 **/
template <typename T>
typename matrix_archetype<T>::size_type num_cols(matrix_archetype<T> const& m) { return typename matrix_archetype<T>::size_type(0); }

} // end namespace concepts
} // end namespace numeric
} // end namespace alps
