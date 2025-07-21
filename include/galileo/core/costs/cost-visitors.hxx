#ifndef __galileo_core_costs_cost_visitors_hxx__
#define __galileo_core_costs_cost_visitors_hxx__

#include <vector>

#include "galileo/core/costs/cost-unary-visitor.hpp"
#include <boost/fusion/container/generation/make_vector.hpp>

#include "galileo/core/costs/cost-visitors.hpp"

#include "galileo/common/container/aligned-vector.hpp"

namespace galileo
{

    // Cost model visitors

    template <typename PhaseSpec, typename StateVectorType, typename ControlVectorType>
    struct CostCalcZerothOrderVisitor
        : fusion::CostUnaryVisitorBase<CostCalcZerothOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType>>
    {
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<ControlVectorType>>;

        template <typename CostModelType>
        static void algo(
            const CostModelBase<CostModelType, PhaseSpec> &cost_model,
            CostDataBase<typename CostModelType::Data_t, PhaseSpec> &cost_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<ControlVectorType> &u)
        {
            cost_model.calc(cost_data, x.derived(), u.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    inline void cost_calc_zeroth_order(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Eigen::MatrixBase<ControlVectorType> &u)
    {
        typedef CostCalcZerothOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType> Algo;

        Algo::run(cost_model, cost_data, typename Algo::ArgsType(x, u));
    }

    template <typename PhaseSpec, typename StateVectorType>
    struct CostCalcZerothOrderVisitor<PhaseSpec, StateVectorType, Blank>
        : fusion::CostUnaryVisitorBase<CostCalcZerothOrderVisitor<PhaseSpec, StateVectorType, Blank>>
    {
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Blank>;

        template <typename CostModelType>
        static void algo(
            const CostModelBase<CostModelType, PhaseSpec> &cost_model,
            CostDataBase<typename CostModelType::Data_t, PhaseSpec> &cost_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Blank blank)
        {
            cost_model.calc(cost_data, x.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename StateVectorType>
    inline void cost_calc_zeroth_order(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Blank blank)
    {
        typedef CostCalcZerothOrderVisitor<PhaseSpec, StateVectorType, Blank> Algo;

        Algo::run(cost_model, cost_data, typename Algo::ArgsType(x, blank));
    }

    template <typename PhaseSpec, typename StateVectorType, typename ControlVectorType>
    struct CostCalcFirstOrderVisitor
        : fusion::CostUnaryVisitorBase<CostCalcFirstOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType>>
    {
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<ControlVectorType>>;

        template <typename CostModelType>
        static void algo(
            const CostModelBase<CostModelType, PhaseSpec> &cost_model,
            CostDataBase<typename CostModelType::Data_t, PhaseSpec> &cost_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<ControlVectorType> &u)
        {
            cost_model.calcDiff(cost_data, x.derived(), u.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    inline void cost_calc_first_order(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Eigen::MatrixBase<ControlVectorType> &u)
    {
        typedef CostCalcFirstOrderVisitor<PhaseSpec, StateVectorType, ControlVectorType> Algo;

        Algo::run(cost_model, cost_data, typename Algo::ArgsType(x, u));
    }

    template <typename PhaseSpec, typename StateVectorType>
    struct CostCalcFirstOrderVisitor<PhaseSpec, StateVectorType, Blank>
        : fusion::CostUnaryVisitorBase<CostCalcFirstOrderVisitor<PhaseSpec, StateVectorType, Blank>>
    {
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Blank>;

        template <typename CostModelType>
        static void algo(
            const CostModelBase<CostModelType, PhaseSpec> &cost_model,
            CostDataBase<typename CostModelType::Data_t, PhaseSpec> &cost_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Blank blank)
        {
            cost_model.calcDiff(cost_data, x.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename StateVectorType>
    inline void cost_calc_first_order(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Blank blank)
    {
        typedef CostCalcFirstOrderVisitor<PhaseSpec, StateVectorType, Blank> Algo;

        Algo::run(cost_model, cost_data, typename Algo::ArgsType(x, blank));
    }

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename DataCollector>
    struct CostCreateDataVisitor
        : fusion::CostUnaryVisitorBase<CostCreateDataVisitor<PhaseSpec, CostCollectionTpl, DataCollector>,
                                       CostDataTpl<PhaseSpec, CostCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<DataCollector *const>;
        using CostCollection_t = CostCollectionTpl<PhaseSpec>;
        using CostModelVariant_t = CostCollection_t::CostModelVariant_t;
        using CostDataVariant_t = CostDataTpl<PhaseSpec, CostCollectionTpl>;

        template <typename CostModelType>
        static CostDataVariant_t algo(
            const CostModelBase<CostModelType, PhaseSpec> &cost_model,
            DataCollector *const collector)
        {
            return CostDataVariant_t(cost_model.createData(collector));
        }
    };

    template <typename PhaseSpec,
              template <typename> class CostCollectionTpl,
              typename DataCollector>
    inline CostDataTpl<PhaseSpec, CostCollectionTpl> cost_create_data(
        const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model,
        DataCollector *const collector)
    {
        typedef CostCreateDataVisitor<PhaseSpec, CostCollectionTpl, DataCollector> Algo;

        return Algo::run(cost_model, typename Algo::ArgsType(collector));
    }

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostGetPhaseSpecVisitor
        : boost::static_visitor<const PhaseSpec &>
    {

        using ReturnType = const PhaseSpec &;

        template <typename CostModelType>
        ReturnType operator()(const CostModelBase<CostModelType, PhaseSpec> &cost_model) const
        {
            return cost_model.get_ps();
        }

        static ReturnType run(const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model)
        {
            return boost::apply_visitor(CostGetPhaseSpecVisitor<PhaseSpec, CostCollectionTpl>(), cost_model);
        }
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    inline const PhaseSpec &cost_get_ps(const CostModelTpl<PhaseSpec, CostCollectionTpl> &cost_model)
    {
        return CostGetPhaseSpecVisitor<PhaseSpec, CostCollectionTpl>::run(cost_model);
    }

    // Cost data visitors

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostLVisitor
        : boost::static_visitor<typename CostDataTpl<PhaseSpec, CostCollectionTpl>::L_t>
    {

        using ReturnType = typename CostDataTpl<PhaseSpec, CostCollectionTpl>::L_t;

        template <typename CostDataType>
        ReturnType operator()(const CostDataBase<CostDataType, PhaseSpec> &cost_data) const
        {
            return cost_data.L();
        }

        static ReturnType run(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
        {
            return boost::apply_visitor(CostLVisitor<PhaseSpec, CostCollectionTpl>(), cost_data);
        }
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::L_t cost_L(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
    {
        return CostLVisitor<PhaseSpec, CostCollectionTpl>::run(cost_data);
    }

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostLxVisitor
        : boost::static_visitor<typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lx_t>
    {
        using ReturnType = typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lx_t;

        template <typename CostDataType>
        ReturnType operator()(const CostDataBase<CostDataType, PhaseSpec> &cost_data) const
        {
            return cost_data.Lx();
        }

        static ReturnType run(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
        {
            return boost::apply_visitor(CostLxVisitor<PhaseSpec, CostCollectionTpl>(), cost_data);
        }
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lx_t cost_Lx(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
    {
        return CostLxVisitor<PhaseSpec, CostCollectionTpl>::run(cost_data);
    }

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostLuVisitor
        : boost::static_visitor<typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lu_t>
    {
        using ReturnType = typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lu_t;

        template <typename CostDataType>
        ReturnType operator()(const CostDataBase<CostDataType, PhaseSpec> &cost_data) const
        {
            return cost_data.Lu();
        }

        static ReturnType run(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
        {
            return boost::apply_visitor(CostLuVisitor<PhaseSpec, CostCollectionTpl>(), cost_data);
        }
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lu_t cost_Lu(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
    {
        return CostLuVisitor<PhaseSpec, CostCollectionTpl>::run(cost_data);
    }

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostLxxVisitor
        : boost::static_visitor<typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lxx_t>
    {
        using ReturnType = typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lxx_t;

        template <typename CostDataType>
        ReturnType operator()(const CostDataBase<CostDataType, PhaseSpec> &cost_data) const
        {
            return cost_data.Lxx();
        }

        static ReturnType run(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
        {
            return boost::apply_visitor(CostLxxVisitor<PhaseSpec, CostCollectionTpl>(), cost_data);
        }
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lxx_t cost_Lxx(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
    {
        return CostLxxVisitor<PhaseSpec, CostCollectionTpl>::run(cost_data);
    }

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostLxuVisitor
        : boost::static_visitor<typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lxu_t>
    {
        using ReturnType = typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lxu_t;

        template <typename CostDataType>
        ReturnType operator()(const CostDataBase<CostDataType, PhaseSpec> &cost_data) const
        {
            return cost_data.Lxu();
        }

        static ReturnType run(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
        {
            return boost::apply_visitor(CostLxuVisitor<PhaseSpec, CostCollectionTpl>(), cost_data);
        }
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Lxu_t cost_Lxu(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
    {
        return CostLxuVisitor<PhaseSpec, CostCollectionTpl>::run(cost_data);
    }

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    struct CostLuuVisitor
        : boost::static_visitor<typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Luu_t>
    {
        using ReturnType = typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Luu_t;

        template <typename CostDataType>
        ReturnType operator()(const CostDataBase<CostDataType, PhaseSpec> &cost_data) const
        {
            return cost_data.Luu();
        }

        static ReturnType run(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
        {
            return boost::apply_visitor(CostLuuVisitor<PhaseSpec, CostCollectionTpl>(), cost_data);
        }
    };

    template <typename PhaseSpec, template <typename> class CostCollectionTpl>
    inline typename CostDataTpl<PhaseSpec, CostCollectionTpl>::Luu_t cost_Luu(const CostDataTpl<PhaseSpec, CostCollectionTpl> &cost_data)
    {
        return CostLuuVisitor<PhaseSpec, CostCollectionTpl>::run(cost_data);
    }

} // namespace galileo

#endif // __galileo_core_costs_cost_visitors_hxx__
