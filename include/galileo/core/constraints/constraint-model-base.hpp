#ifndef __galileo_core_constraints_constraint_model_base_hpp__
#define __galileo_core_constraints_constraint_model_base_hpp__

#include "galileo/core/constraints/constraint-base.hpp"

// #define GALILEO_CONSTRAINT_BASIC_TYPEDEF(Constraint)                                    \
//     using VarScalar = typename traits<Constraint>::VarScalar;                           \
//     using NumScalar = typename traits<Constraint>::NumScalar;                           \
//     static constexpr int Options = traits<Constraint>::Options;                         \
//     using ConstraintModelDerived = typename traits<Constraint>::ConstraintModelDerived; \
//     using ConstraintDataDerived = typename traits<Constraint>::ConstraintDataDerived;

// #define GALILEO_CONSTRAINT_CONSTANTS(Constraint)      \
//     static constexpr int NX = traits<Constraint>::NX; \
//     static constexpr int NU = traits<Constraint>::NU; \
//     static constexpr int NH = traits<Constraint>::NH; \
//     static constexpr int NG = traits<Constraint>::NG;

// #define GALILEO_CONSTRAINT_MODEL_TYPEDEF(Constraint) \
//     using ResidualModel_t = typename traits<Constraint>::ResidualModel_t;

// #define GALILEO_CONSTRAINT_DATA_TYPEDEF(Constraint)                     \
//     using ResidualData_t = typename traits<Constraint>::ResidualData_t; \
//     using H_t = typename traits<Constraint>::H_t;                       \
//     using Hx_t = typename traits<Constraint>::Hx_t;                     \
//     using Hu_t = typename traits<Constraint>::Hu_t;                     \
//     using G_t = typename traits<Constraint>::G_t;                       \
//     using Gx_t = typename traits<Constraint>::Gx_t;                     \
//     using Gu_t = typename traits<Constraint>::Gu_t;

namespace galileo
{
    namespace core
    {

        template <typename Derived, typename PhaseSpec>
        class ConstraintModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using ConstraintDerived = typename traits<Derived>::ConstraintDerived;
            using ConstraintModelDerived = typename traits<ConstraintDerived>::ConstraintModelDerived;
            using ConstraintDataDerived = typename traits<ConstraintDerived>::ConstraintDataDerived;

            using BoundVector_t = Eigen::Matrix<typename PS::NumScalar, traits<Derived>::NG, 1, PS::Options>;

            template <typename StateVectorType, typename ControlVectorType>
            void calc(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calc(data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calc(ConstraintDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calc(data, x.derived());
            }

            template <typename StateVectorType, typename ControlVectorType>
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
            {
                derived().calcDiff(data, x.derived(), u.derived());
            }

            template <typename StateVectorType>
            void calcDiff(ConstraintDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
            {
                derived().calcDiff(data, x.derived());
            }

            template <typename DataCollector>
            auto createData(DataCollector *const collector)
            {
                return derived().createData(collector);
            }

            template <typename LowerBoundType, typename UpperBoundType>
            void updateBounds(const Eigen::MatrixBase<LowerBoundType> &lb,
                              const Eigen::MatrixBase<UpperBoundType> &ub)
            {
                derived().updateBounds(lb.derived(), ub.derived());
            }

            void removeBounds()
            {
                derived().removeBounds();
            }

            const BoundVector_t &lb() const
            {
                return derived().lb_impl();
            }

            const BoundVector_t &ub() const
            {
                return derived().ub_impl();
            }

            int ng() const
            {
                return derived().ng_impl();
            }

            int nh() const
            {
                return derived().nh_impl();
            }

            int ng_impl() const
            {
                return traits<Derived>::NG;
            }

            int nh_impl() const
            {
                return traits<Derived>::NH;
            }

        protected:
            inline ConstraintModelBase()
            {
            }

            inline ConstraintModelBase(const ConstraintModelBase &clone)
            {
                *this = clone;
            }

            inline ConstraintModelBase &operator=(const ConstraintModelBase &clone)
            {
                return *this;
            }

        }; // class ConstraintModelBase

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_model_base_hpp__
