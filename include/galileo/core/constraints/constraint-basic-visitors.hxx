#ifndef __galileo_core_constraints_constraint_basic_visitors_hxx__
#define __galileo_core_constraints_constraint_basic_visitors_hxx__

#include <vector>

#include <boost/fusion/container/generation/make_vector.hpp>
#include "galileo/core/constraints/constraint-unary-visitor.hpp"

#include "galileo/core/constraints/constraint-basic-visitors.hpp"

#include "galileo/utils/aligned-vector.hpp"

namespace galileo
{

    namespace core
    {

        // Constraint model visitors

        template <typename PhaseSpec, typename StateVectorType, typename ControlVectorType>
        struct ConstraintCalcZerothOrderVisitor
            : fusion::ConstraintUnaryVisitorBase<ConstraintCalcZerothOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<ControlVectorType>>;

            template <typename ConstraintModel>
            static void algo(
                const core::ConstraintModelBase<ConstraintModel, PhaseSpec> &constraint_model,
                typename core::ConstraintDataBase<typename ConstraintModel::ConstraintDataDerived, PhaseSpec> &constraint_data,
                const Eigen::MatrixBase<StateVectorType> &x,
                const Eigen::MatrixBase<ControlVectorType> &u)
            {
                constraint_model.calc(constraint_data, x.derived(), u.derived());
            }
        };

        template <typename PhaseSpec,
                  template <typename> class ConstraintCollectionTpl,
                  typename StateVectorType,
                  typename ControlVectorType>
        inline void constraint_calc_zeroth_order(
            const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
            ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<ControlVectorType> &u)
        {
            typedef ConstraintCalcZerothOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType> Algo;

            Algo::run(constraint_model, constraint_data, typename Algo::ArgsType(x, u));
        }

        template <typename PhaseSpec, typename StateVectorType, typename ControlVectorType>
        struct ConstraintCalcFirstOrderVisitor
            : fusion::ConstraintUnaryVisitorBase<ConstraintCalcFirstOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<ControlVectorType>>;

            template <typename ConstraintModel>
            static void algo(
                const core::ConstraintModelBase<ConstraintModel, PhaseSpec> &constraint_model,
                typename core::ConstraintDataBase<typename ConstraintModel::ConstraintDataDerived, PhaseSpec> &constraint_data,
                const Eigen::MatrixBase<StateVectorType> &x,
                const Eigen::MatrixBase<ControlVectorType> &u)
            {
                constraint_model.calcDiff(constraint_data, x.derived(), u.derived());
            }
        };

        template <typename PhaseSpec,
                  template <typename> class ConstraintCollectionTpl,
                  typename StateVectorType,
                  typename ControlVectorType>
        inline void constraint_calc_first_order(
            const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
            ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<ControlVectorType> &u)
        {
            typedef ConstraintCalcFirstOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType> Algo;

            Algo::run(constraint_model, constraint_data, typename Algo::ArgsType(x, u));
        }

        template <typename PhaseSpec, typename LowerBoundType, typename UpperBoundType>
        struct ConstraintUpdateBoundsVisitor
            : fusion::ConstraintUnaryVisitorBase<ConstraintUpdateBoundsVisitor<PhaseSpec, LowerBoundType, UpperBoundType>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<LowerBoundType>, Eigen::MatrixBase<UpperBoundType>>;

            template <typename ConstraintModel>
            static void algo(
                const core::ConstraintModelBase<ConstraintModel, PhaseSpec> &constraint_model,
                const Eigen::MatrixBase<LowerBoundType> &lb,
                const Eigen::MatrixBase<UpperBoundType> &ub)
            {
                constraint_model.updateBounds(lb.derived(), ub.derived());
            }
        };

        template <typename PhaseSpec,
                  template <typename> class ConstraintCollectionTpl,
                  typename LowerBoundType,
                  typename UpperBoundType>
        inline void constraint_update_bounds(
            const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
            const Eigen::MatrixBase<LowerBoundType> &lb,
            const Eigen::MatrixBase<UpperBoundType> &ub)
        {
            typedef ConstraintUpdateBoundsVisitor<PhaseSpec, LowerBoundType, UpperBoundType> Algo;

            Algo::run(constraint_model, typename Algo::ArgsType(lb, ub));
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintLbVisitor : boost::static_visitor<typename ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>::BoundVector_t>
        {
            template <typename ConstraintModelDerived>
            typename ConstraintModelDerived::BoundVector_t operator()(const ConstraintModelBase<ConstraintModelDerived, PhaseSpec> &constraint_model) const
            {
                return constraint_model.lb();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static typename ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>::BoundVector_t run(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
            {
                return boost::apply_visitor(ConstraintLbVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_model);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline typename ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>::BoundVector_t constraint_lb(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
        {
            return ConstraintLbVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_model);
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintUbVisitor : boost::static_visitor<typename ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>::BoundVector_t>
        {
            template <typename ConstraintModelDerived>
            typename ConstraintModelDerived::BoundVector_t operator()(const ConstraintModelBase<ConstraintModelDerived, PhaseSpec> &constraint_model) const
            {
                return constraint_model.ub();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static typename ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>::BoundVector_t run(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
            {
                return boost::apply_visitor(ConstraintUbVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_model);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline typename ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl>::BoundVector_t constraint_ub(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
        {
            return ConstraintUbVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_model);
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintNgVisitor : boost::static_visitor<int>
        {
            template <typename ConstraintModelDerived>
            int operator()(const ConstraintModelBase<ConstraintModelDerived, PhaseSpec> &constraint_model) const
            {
                return constraint_model.ng();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static int run(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
            {
                return boost::apply_visitor(ConstraintNgVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_model);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline int constraint_ng(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
        {
            return ConstraintNgVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_model);
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintNhVisitor : boost::static_visitor<int>
        {
            template <typename ConstraintModelDerived>
            int operator()(const ConstraintModelBase<ConstraintModelDerived, PhaseSpec> &constraint_model) const
            {
                return constraint_model.nh();
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline int constraint_nh(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
        {
            return ConstraintNhVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_model);
        }

        // Constraint data visitors

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintHVisitor : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::H_t>
        {

            using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::H_t;

            template <typename ConstraintDataDerived>
            ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived, PhaseSpec> &constraint_data) const
            {
                return constraint_data.H();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static ReturnType run(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
            {
                return boost::apply_visitor(ConstraintHVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_data);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::H_t constraint_H(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
        {
            return ConstraintHVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_data);
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintHxVisitor : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hx_t>
        {
            using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hx_t;

            template <typename ConstraintDataDerived>
            ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived, PhaseSpec> &constraint_data) const
            {
                return constraint_data.Hx();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static ReturnType run(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
            {
                return boost::apply_visitor(ConstraintHxVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_data);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hx_t constraint_Hx(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
        {
            return ConstraintHxVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_data);
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintHuVisitor : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hu_t>
        {
            using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hu_t;

            template <typename ConstraintDataDerived>
            ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived, PhaseSpec> &constraint_data) const
            {
                return constraint_data.Hu();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static ReturnType run(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
            {
                return boost::apply_visitor(ConstraintHuVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_data);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hu_t constraint_Hu(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
        {
            return ConstraintHuVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_data);
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintGVisitor : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::G_t>
        {
            using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::G_t;

            template <typename ConstraintDataDerived>
            ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived, PhaseSpec> &constraint_data) const
            {
                return constraint_data.G();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static ReturnType run(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
            {
                return boost::apply_visitor(ConstraintGVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_data);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::G_t constraint_G(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
        {
            return ConstraintGVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_data);
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintGxVisitor : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Gx_t>
        {
            using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Gx_t;

            template <typename ConstraintDataDerived>
            ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived, PhaseSpec> &constraint_data) const
            {
                return constraint_data.Gx();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static ReturnType run(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
            {
                return boost::apply_visitor(ConstraintGxVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_data);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Gx_t constraint_Gx(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
        {
            return ConstraintGxVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_data);
        }

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        struct ConstraintGuVisitor : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Gu_t>
        {
            using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Gu_t;

            template <typename ConstraintDataDerived>
            ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived, PhaseSpec> &constraint_data) const
            {
                return constraint_data.Gu();
            }

            template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
            static ReturnType run(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
            {
                return boost::apply_visitor(ConstraintGuVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_data);
            }
        };

        template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
        inline typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Gu_t constraint_Gu(const ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data)
        {
            return ConstraintGuVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_data);
        }

    } // namespace core

} // namespace galileo

#endif // __galileo_core_constraints_constraint_basic_visitors_hxx__