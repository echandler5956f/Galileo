#ifndef __galileo_predictive_phases_phase_visitors_hxx__
#define __galileo_predictive_phases_phase_visitors_hxx__

#include <vector>

#include "galileo/predictive/phases/phase-unary-visitor.hpp"
#include <boost/fusion/container/generation/make_vector.hpp>

#include "galileo/predictive/phases/phase-visitors.hpp"

#include "galileo/common/container/aligned-vector.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename> class PhaseCollectionTpl,
              typename DataCollector>
    struct PhaseCreateDataVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCreateDataVisitor<PhaseSpec, PhaseCollectionTpl, DataCollector>>
    {
        using ArgsType = boost::fusion::vector<DataCollector *const>;
        using PhaseCollection_t = PhaseCollectionTpl<PhaseSpec>;
        using PhaseModelVariant_t = PhaseCollection_t::ModelVariant_t;
        using PhaseDataVariant_t = PhaseDataTpl<PhaseSpec, PhaseCollectionTpl>;

        template <typename PhaseModelDerived>
        static PhaseDataVariant_t algo(
            const PhaseModelBase<PhaseModelDerived, PhaseSpec> &phase_model,
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

    template <typename StateMatrixType, typename ControlParamMatrixType>
    struct PhaseCalcZerothOrderVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCalcZerothOrderVisitor<StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateMatrixType>, Eigen::MatrixBase<ControlParamMatrixType>>;

        template <typename PhaseModel>
        static void algo(
            const PhaseModelBase<PhaseModel> &phase_model,
            PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            const Eigen::MatrixBase<ControlParamMatrixType> &ws)
        {
            phase_model.calc(phase_data, xs.derived(), ws.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_zeroth_order(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws)
    {
        typedef PhaseCalcZerothOrderVisitor<StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs, ws));
    }

    template <typename StateMatrixType, typename ControlParamMatrixType>
    struct PhaseCalcFirstOrderVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCalcFirstOrderVisitor<StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateMatrixType>, Eigen::MatrixBase<ControlParamMatrixType>>;

        template <typename PhaseModel>
        static void algo(
            const galileo::PhaseModelBase<PhaseModel> &phase_model,
            typename galileo::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            const Eigen::MatrixBase<ControlParamMatrixType> &ws)
        {
            phase_model.calcDiff(phase_data, xs.derived(), ws.derived());
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_first_order(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws)
    {
        typedef PhaseCalcFirstOrderVisitor<StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs, ws));
    }

    template <
        typename PhaseSpec,
        template <typename PS> class PhaseCollectionTpl,
        typename StateMatrixType,
        typename ControlParamMatrixType>
    struct PhaseQuasiStaticVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseQuasiStaticVisitor<PhaseSpec, PhaseCollectionTpl, StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateMatrixType>, Eigen::MatrixBase<ControlParamMatrixType>, const std::size_t &, const typename PhaseSpec::NumScalar &>;

        template <typename PhaseModel>
        static void algo(
            const galileo::PhaseModelBase<PhaseModel> &phase_model,
            typename galileo::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            Eigen::MatrixBase<ControlParamMatrixType> &ws,
            const std::size_t &maxiter,
            const typename PhaseSpec::NumScalar &tol)
        {
            phase_model.quasiStatic(phase_data, xs.derived(), ws.derived(), maxiter, tol);
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_quasi_static(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        Eigen::MatrixBase<ControlParamMatrixType> &ws,
        const std::size_t &maxiter,
        const typename PhaseSpec::NumScalar &tol)
    {
        typedef PhaseQuasiStaticVisitor<PhaseSpec, PhaseCollectionTpl, StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs, ws, maxiter, tol));
    }

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    struct PhaseSegmentModelVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseSegmentModelVisitor<PhaseSpec, PhaseCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<>;

        template <typename PhaseModel>
        static typename PhaseSpec::SegmentModel_t algo(
            const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model)
        {
            return phase_model.segment();
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    inline typename PhaseSpec::SegmentModel_t phase_segment_model(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model)
    {
        typedef PhaseSegmentModelVisitor<PhaseSpec, PhaseCollectionTpl> Algo;

        return Algo::run(phase_model, typename Algo::ArgsType());
    }

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    struct PhasePeriodVisitor
        : fusion::PhaseUnaryVisitorBase<PhasePeriodVisitor<PhaseSpec, PhaseCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<>;

        template <typename PhaseModel>
        static typename PhaseSpec::NumScalar algo(
            const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model)
        {
            return phase_model.period();
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    inline typename PhaseSpec::NumScalar phase_period(
        const PhaseModelTpl<PhaseSpec, PhaseCollectionTpl> &phase_model)
    {
        typedef PhasePeriodVisitor<PhaseSpec, PhaseCollectionTpl> Algo;

        return Algo::run(phase_model, typename Algo::ArgsType());
    }

    // Phase data visitors

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    struct PhaseSegmentDataVectorVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseSegmentDataVectorVisitor<PhaseSpec, PhaseCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<>;

        template <typename PhaseData>
        static typename PhaseSpec::SegmentDataVector_t algo(
            PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data)
        {
            return phase_data.segments();
        }
    };

    template <typename PhaseSpec,
              template <typename PS> class PhaseCollectionTpl>
    inline typename PhaseSpec::SegmentDataVector_t phase_segment_data_vector(
        PhaseDataTpl<PhaseSpec, PhaseCollectionTpl> &phase_data)
    {
        typedef PhaseSegmentDataVectorVisitor<PhaseSpec, PhaseCollectionTpl> Algo;

        return Algo::run(phase_data, typename Algo::ArgsType());
    }

} // namespace galileo

#endif // __galileo_predictive_phases_phase_visitors_hxx__
