#ifndef __galileo_utils_aligned_vector_hpp__
#define __galileo_utils_aligned_vector_hpp__

#include <vector>
#include <Eigen/StdVector>

#define GALILEO_ALIGNED_STD_VECTOR(Type) ::galileo::container::aligned_vector<Type>
#define GALILEO_STD_VECTOR_WITH_EIGEN_ALLOCATOR(T) ::galileo::container::aligned_vector<T>

namespace galileo
{
    namespace container
    {

        template <typename T>
        using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;

    } // namespace container

} // namespace galileo

#endif // __galileo_utils_aligned_vector_hpp__