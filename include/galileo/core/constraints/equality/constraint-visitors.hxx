#ifndef __galileo_core_constraints_equality_constraint_visitors_hxx__
#define __galileo_core_constraints_equality_constraint_visitors_hxx__

#include <vector>

#include "galileo/core/constraints/equality/constraint-unary-visitor.hpp"
#include <boost/fusion/container/generation/make_vector.hpp>

#include "galileo/core/constraints/equality/constraint-visitors.hpp"

#include "galileo/common/container/aligned-vector.hpp"

namespace galileo
{

    // Constraint model visitors

    template <typename PhaseSpec, typename StateVectorType, typename ControlVectorType>
    struct ConstraintCalcZerothOrderVisitor
        : fusion::ConstraintUnaryVisitorBase<ConstraintCalcZerothOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &, const ControlVectorType &>;

        template <typename ConstraintModelType>
        static void algo(
            const ConstraintModelBase<ConstraintModelType, PhaseSpec> &constraint_model,
            ConstraintDataBase<typename ConstraintModelType::Data_t, PhaseSpec> &constraint_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<ControlVectorType> &u)
        {
            constraint_model.calc(constraint_data.derived(), x.derived(), u.derived());
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

        Algo::run(constraint_model, constraint_data, typename Algo::ArgsType(x.derived(), u.derived()));
    }

    template <typename PhaseSpec, typename StateVectorType>
    struct ConstraintCalcZerothOrderVisitor<PhaseSpec, StateVectorType, Blank>
        : fusion::ConstraintUnaryVisitorBase<ConstraintCalcZerothOrderVisitor<PhaseSpec, StateVectorType, Blank>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &, Blank>;

        template <typename ConstraintModelType>
        static void algo(
            const ConstraintModelBase<ConstraintModelType, PhaseSpec> &constraint_model,
            ConstraintDataBase<typename ConstraintModelType::Data_t, PhaseSpec> &constraint_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Blank blank)
        {
            constraint_model.calc(constraint_data.derived(), x.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename StateVectorType>
    inline void constraint_calc_zeroth_order(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Blank blank)
    {
        typedef ConstraintCalcZerothOrderVisitor<PhaseSpec, StateVectorType, Blank> Algo;

        Algo::run(constraint_model, constraint_data, typename Algo::ArgsType(x.derived(), blank));
    }

    template <typename PhaseSpec, typename StateVectorType, typename ControlVectorType>
    struct ConstraintCalcFirstOrderVisitor
        : fusion::ConstraintUnaryVisitorBase<ConstraintCalcFirstOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &, const ControlVectorType &>;

        template <typename ConstraintModelType>
        static void algo(
            const ConstraintModelBase<ConstraintModelType, PhaseSpec> &constraint_model,
            ConstraintDataBase<typename ConstraintModelType::Data_t, PhaseSpec> &constraint_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<ControlVectorType> &u)
        {
            constraint_model.calcDiff(constraint_data.derived(), x.derived(), u.derived());
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

        Algo::run(constraint_model, constraint_data, typename Algo::ArgsType(x.derived(), u.derived()));
    }

    template <typename PhaseSpec, typename StateVectorType>
    struct ConstraintCalcFirstOrderVisitor<PhaseSpec, StateVectorType, Blank>
        : fusion::ConstraintUnaryVisitorBase<ConstraintCalcFirstOrderVisitor<PhaseSpec, StateVectorType, Blank>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &, Blank>;

        template <typename ConstraintModelType>
        static void algo(
            const ConstraintModelBase<ConstraintModelType, PhaseSpec> &constraint_model,
            ConstraintDataBase<typename ConstraintModelType::Data_t, PhaseSpec> &constraint_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Blank blank)
        {
            constraint_model.calcDiff(constraint_data.derived(), x.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename StateVectorType>
    inline void constraint_calc_first_order(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Blank blank)
    {
        typedef ConstraintCalcFirstOrderVisitor<PhaseSpec, StateVectorType, Blank> Algo;

        Algo::run(constraint_model, constraint_data, typename Algo::ArgsType(x.derived(), blank));
    }

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename DataCollector>
    struct ConstraintCreateDataVisitor
        : fusion::ConstraintUnaryVisitorBase<ConstraintCreateDataVisitor<PhaseSpec, ConstraintCollectionTpl, DataCollector>,
                                             ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<DataCollector *const>;
        using ConstraintCollection_t = ConstraintCollectionTpl<PhaseSpec>;
        using ConstraintModelVariant_t = ConstraintCollection_t::ConstraintModelVariant_t;
        using ConstraintDataVariant_t = ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>;

        template <typename ConstraintModelType>
        static ConstraintDataVariant_t algo(
            const ConstraintModelBase<ConstraintModelType, PhaseSpec> &constraint_model,
            DataCollector *const collector)
        {
            return ConstraintDataVariant_t(constraint_model.createData(collector));
        }
    };

    template <typename PhaseSpec,
              template <typename> class ConstraintCollectionTpl,
              typename DataCollector>
    inline ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl> constraint_create_data(
        const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model,
        DataCollector *const collector)
    {
        typedef ConstraintCreateDataVisitor<PhaseSpec, ConstraintCollectionTpl, DataCollector> Algo;

        return Algo::run(constraint_model, typename Algo::ArgsType(collector));
    }

    template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
    struct ConstraintGetPhaseSpecVisitor
        : boost::static_visitor<const PhaseSpec &>
    {

        using ReturnType = PhaseSpec;

        template <typename ConstraintModelType>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelType, PhaseSpec> &constraint_model) const
        {
            return constraint_model.get_ps();
        }

        static ReturnType run(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
        {
            return boost::apply_visitor(ConstraintGetPhaseSpecVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_model);
        }
    };

    template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
    inline PhaseSpec constraint_get_ps(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
    {
        return ConstraintGetPhaseSpecVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_model);
    }

    template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
    struct ConstraintGetNHVisitor
        : boost::static_visitor<int>
    {

        using ReturnType = int;

        template <typename ConstraintModelType>
        ReturnType operator()(const ConstraintModelBase<ConstraintModelType, PhaseSpec> &constraint_model) const
        {
            return constraint_model.get_nh();
        }

        ReturnType operator()(const ConstraintModelVoid &) const
        {
            return 0;
        }

        static ReturnType run(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
        {
            return boost::apply_visitor(ConstraintGetNHVisitor<PhaseSpec, ConstraintCollectionTpl>(), constraint_model);
        }
    };

    template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
    inline int constraint_get_nh(const ConstraintModelTpl<PhaseSpec, ConstraintCollectionTpl> &constraint_model)
    {
        return ConstraintGetNHVisitor<PhaseSpec, ConstraintCollectionTpl>::run(constraint_model);
    }

    // Constraint data visitors

    template <typename PhaseSpec, template <typename> class ConstraintCollectionTpl>
    struct ConstraintHVisitor
        : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::H_t>
    {

        using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::H_t;

        template <typename ConstraintDataType>
        ReturnType operator()(const ConstraintDataBase<ConstraintDataType, PhaseSpec> &constraint_data) const
        {
            return constraint_data.H();
        }

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
    struct ConstraintHxVisitor
        : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hx_t>
    {
        using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hx_t;

        template <typename ConstraintDataDerived>
        ReturnType operator()(const ConstraintDataBase<ConstraintDataDerived, PhaseSpec> &constraint_data) const
        {
            return constraint_data.Hx();
        }

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
    struct ConstraintHuVisitor
        : boost::static_visitor<typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hu_t>
    {
        using ReturnType = typename ConstraintDataTpl<PhaseSpec, ConstraintCollectionTpl>::Hu_t;

        template <typename ConstraintDataType>
        ReturnType operator()(const ConstraintDataBase<ConstraintDataType, PhaseSpec> &constraint_data) const
        {
            return constraint_data.Hu();
        }

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

} // namespace galileo

#endif // __galileo_core_constraints_equality_constraint_visitors_hxx__
