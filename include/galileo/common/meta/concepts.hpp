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

    enum AssignmentOp
    {
        SETTO,
        ADDTO,
        RMFROM
    }; // enum AssignmentOp

    template <AssignmentOp op>
    concept IsSetTo = (op == SETTO);
    template <AssignmentOp op>
    concept IsAddTo = (op == ADDTO);
    template <AssignmentOp op>
    concept IsRmFrom = (op == RMFROM);

    // generic, but intentionally lines up with pinocchio::ArgumentPosition
    enum Jcomponent
    {
        FIRST = 0,
        SECOND = 1,
        BOTH = 2
    }; // enum Jcomponent

    template <Jcomponent jc>
    concept IsFirst = (jc == FIRST);
    template <Jcomponent jc>
    concept IsSecond = (jc == SECOND);
    template <Jcomponent jc>
    concept IsBoth = (jc == BOTH);

} // namespace galileo

#endif // __galileo_common_meta_concepts_hpp__
