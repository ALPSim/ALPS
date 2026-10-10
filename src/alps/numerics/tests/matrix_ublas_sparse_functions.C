#include "matrix_unit_tests.hpp"
#include <alps/numeric/matrix/ublas_sparse_functions.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

template<class T> class SparseMatrix : public ::testing::Test {};
TYPED_TEST_SUITE(SparseMatrix, test_types);

TYPED_TEST(SparseMatrix, sparse_matrix_vector_multiply)
{
    using T = TypeParam;
    namespace ublas = boost::numeric::ublas;
    alps::numeric::vector<T> v(10);
    fill_range_with_numbers(v.begin(),v.end(),0);

    ublas::mapped_vector_of_mapped_vector<T,ublas::row_major> m(20,10);

    for(int i=0; i < 9; ++i)
    {
        m.insert_element(2*i,i,10*i);
        m.insert_element(2*i,i+1,100*i);
    }

    alps::numeric::vector<T> r = m * v;

    std::cout << r << std::endl;
    for(int i=0; i < 9; ++i)
    {
        EXPECT_EQ(r(2*i+1), T(0));
        EXPECT_EQ(r(2*i),T(10*i)*v(i)+T(100*i)*v(i+1));
    }
}
