#ifndef __galileo_utils_cast_hpp_
#define __galileo_utils_cast_hpp_

#include <Eigen/Core>

namespace galileo
{
    
    template <typename NewScalar, typename Scalar>
    NewScalar cast(const Scalar &value)
    {
        return Eigen::internal::cast_impl<Scalar, NewScalar>::run(value);
    }

} // namespace galileo

#endif // __galileo_utils_cast_hpp_