#ifndef __galileo_common_meta_concepts_hpp__
#define __galileo_common_meta_concepts_hpp__

#include <concepts>

namespace galileo
{
    // Concepts to constrain templates to Eigen types
    template <typename Derived>
    concept IsEigenRowVector = (Derived::RowsAtCompileTime == 1 && Derived::ColsAtCompileTime != 1);

    template <typename Derived>
    concept IsEigenColVector = (Derived::ColsAtCompileTime == 1 && Derived::RowsAtCompileTime != 1);

    template <typename Derived>
    concept IsEigenVector = IsEigenRowVector<Derived> || IsEigenColVector<Derived>;

    template <typename Derived>
    concept IsEigenMatrix = !IsEigenVector<Derived>;

} // namespace galileo

#endif // __galileo_common_meta_concepts_hpp__
