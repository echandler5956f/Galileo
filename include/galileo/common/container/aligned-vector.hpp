#ifndef __galileo_common_container_aligned_vector_hpp__
#define __galileo_common_container_aligned_vector_hpp__

#include <Eigen/StdVector>
#include <vector>

#define GALILEO_ALIGNED_STD_VECTOR(Type) ::galileo::aligned_vector<Type>

namespace galileo
{
    template <typename T>
    using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;

} // namespace galileo

#endif // __galileo_common_container_aligned_vector_hpp__
