#ifndef __galileo_common_meta_concepts_hpp__
#define __galileo_common_meta_concepts_hpp__

#include <concepts>

namespace galileo
{
    // Concepts to constrain templates to Eigen types
    template <typename T>
    concept IsEigenDenseBase = std::is_base_of_v<Eigen::DenseBase<std::decay_t<T>>, std::decay_t<T>>;

    template <typename T>
    concept IsEigenMatrixBase = std::is_base_of_v<Eigen::MatrixBase<std::decay_t<T>>, std::decay_t<T>>;

    template <typename T>
    concept IsEigenRowVector = IsEigenMatrixBase<T> && (T::RowsAtCompileTime == 1 && T::ColsAtCompileTime != 1);

    template <typename T>
    concept IsEigenColVector = IsEigenMatrixBase<T> && (T::ColsAtCompileTime == 1 && T::RowsAtCompileTime != 1);

    template <typename T>
    concept IsEigenVector = IsEigenRowVector<T> || IsEigenColVector<T>;

    template <typename T>
    concept IsEigenMatrix = !IsEigenVector<T>;

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
