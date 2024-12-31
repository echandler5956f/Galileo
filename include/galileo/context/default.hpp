#ifndef __galileo_context_default_hpp__
#define __galileo_context_default_hpp__

#define GALILEO_SCALAR_TYPE double

// #include "galileo/context/generic.hpp"

#include <Eigen/Core>
#include "galileo/utils/aligned-vector.hpp"

namespace galileo
{

    template <typename _Scalar, int _Options>
    struct JointCollectionDefaultTpl;

    namespace context
    {
        typedef GALILEO_SCALAR_TYPE Scalar;
        enum
        {
            Options = 0
        };

    } // namespace context

    // Read and write
    template <typename Derived>
    Eigen::Ref<typename Derived::PlainObject> make_ref(const Eigen::MatrixBase<Derived> &x)
    {
        return Eigen::Ref<typename Derived::PlainObject>(x.const_cast_derived());
    }

    // Read-only
    template <typename M>
    auto make_const_ref(Eigen::MatrixBase<M> const &m) -> Eigen::Ref<typename M::PlainObject const>
    {
        return m;
    }

} // namespace galileo

#undef GALILEO_SCALAR_TYPE
#endif // __galileo_context_default_hpp__
