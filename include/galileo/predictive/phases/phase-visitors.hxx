#ifndef __galileo_predictive_phases_phase_visitors_hxx__
#define __galileo_predictive_phases_phase_visitors_hxx__

#include "galileo/predictive/phases/phase-visitor-base.hpp"
#include "galileo/predictive/phases/phase-visitors.hpp"


namespace galileo
{

    template <typename PhaseSpec, typename StateMatrixType, typename ControlParamMatrixType>
    struct PhaseCalcZerothOrderVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCalcZerothOrderVisitor<PhaseSpec, StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<const StateMatrixType &, const ControlParamMatrixType &>;

        template <typename PhaseModelType>
        static void algo(
            const PhaseModelBase<PhaseModelType, PhaseSpec> &phase_model,
            PhaseDataBase<typename PhaseModelType::Data_t, PhaseSpec> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            const Eigen::MatrixBase<ControlParamMatrixType> &ws)
        {
            phase_model.calc(phase_data.derived(), xs.derived(), ws.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_zeroth_order(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws)
    {
        typedef PhaseCalcZerothOrderVisitor<PhaseSpec, StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs.derived(), ws.derived()));
    }

    template <typename PhaseSpec, typename StateMatrixType, typename ControlParamMatrixType>
    struct PhaseCalcFirstOrderVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCalcFirstOrderVisitor<PhaseSpec, StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<const StateMatrixType &, const ControlParamMatrixType &>;

        template <typename PhaseModelType>
        static void algo(
            const PhaseModelBase<PhaseModelType, PhaseSpec> &phase_model,
            PhaseDataBase<typename PhaseModelType::Data_t, PhaseSpec> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            const Eigen::MatrixBase<ControlParamMatrixType> &ws)
        {
            phase_model.calcDiff(phase_data.derived(), xs.derived(), ws.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_first_order(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws)
    {
        typedef PhaseCalcFirstOrderVisitor<PhaseSpec, StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs.derived(), ws.derived()));
    }

    template <typename PhaseSpec,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    struct PhaseQuasiStaticVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseQuasiStaticVisitor<PhaseSpec, StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<const StateMatrixType &, ControlParamMatrixType &, const int, const typename PhaseSpec::NumScalar>;

        template <typename PhaseModelType>
        static void algo(
            const PhaseModelBase<PhaseModelType, PhaseSpec> &phase_model,
            PhaseDataBase<typename PhaseModelType::Data_t, PhaseSpec> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            Eigen::MatrixBase<ControlParamMatrixType> &ws,
            const int maxiter,
            const typename PhaseSpec::NumScalar tol)
        {
            phase_model.quasiStatic(phase_data.derived(), xs.derived(), ws.derived(), maxiter, tol);
        }
    };

    template <typename PhaseSpec,
              template <typename> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_quasi_static(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        Eigen::MatrixBase<ControlParamMatrixType> &ws,
        const int maxiter,
        const typename PhaseSpec::NumScalar tol)
    {
        typedef PhaseQuasiStaticVisitor<PhaseSpec, StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs.derived(), ws.derived(), maxiter, tol));
    }

    template <typename PhaseSpec,
              template <typename> class PhaseCollectionTpl,
              typename DataCollector>
    struct PhaseCreateDataVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCreateDataVisitor<PhaseSpec, PhaseCollectionTpl, DataCollector>,
                                        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<DataCollector *const>;
        using PhaseCollection_t = PhaseCollectionTpl<PhaseSpec>;
        using PhaseModelVariant_t = PhaseCollection_t::PhaseModelVariant_t;
        using PhaseDataVariant_t = PhaseDataTpl<PhaseSpec, PhaseCollectionTpl>;

        template <typename PhaseModelType>
        static PhaseDataVariant_t algo(
            const PhaseModelBase<PhaseModelType, PhaseSpec> &phase_model,
            DataCollector *const collector)
        {
            return PhaseDataVariant_t(phase_model.createData(collector));
        }
    };

    template <typename PhaseSpec,
              template <typename> class PhaseCollectionTpl,
              typename DataCollector>
    inline PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> phase_create_data(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        DataCollector *const collector)
    {
        typedef PhaseCreateDataVisitor<PhaseSpec, PhaseCollectionTpl, DataCollector> Algo;

        return Algo::run(phase_model, typename Algo::ArgsType(collector));
    }

} // namespace galileo

#endif // __galileo_predictive_phases_phase_visitors_hxx__
