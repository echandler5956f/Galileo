#ifndef __galileo_predictive_phases_phase_basic_visitors_hxx__
#define __galileo_predictive_phases_phase_basic_visitors_hxx__

#include <vector>

#include <boost/fusion/container/generation/make_vector.hpp>
#include "galileo/predictive/visitor/phase-unary-visitor.hpp"

#include "galileo/predictive/phases/phase-basic-visitors.hpp"

#include "galileo/utils/aligned-vector.hpp"

namespace galileo
{

    namespace predictive
    {

        template <typename StateVectorType, typename ControlVectorType>
        struct SegmentCalcZerothOrderVisitor
            : fusion::PhaseUnaryVisitorBase<SegmentCalcZerothOrderVisitor<StateVectorType, ControlVectorType>>
        {
            using ArgsType = boost::fusion::vector<std::size_t, Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<ControlVectorType>>;

            template <typename PhaseModel>
            static void algo(
                const predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const std::size_t &segment_index,
                const Eigen::MatrixBase<StateVectorType> &xs,
                const Eigen::MatrixBase<ControlVectorType> &us)
            {
                phase_model.segments_[segment_index].calc(phase_data[segment_index], xs.derived(), us.derived());
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename ControlVectorType>
        inline void segment_calc_zeroth_order(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &segment_index,
            const Eigen::MatrixBase<StateVectorType> &xs,
            const Eigen::MatrixBase<ControlVectorType> &us)
        {
            typedef SegmentCalcZerothOrderVisitor<StateVectorType, ControlVectorType> Algo;

            Algo::run(phase_model, phase_data, typename Algo::ArgsType(segment_index, xs, us));
        }

        template <typename StateVectorType, typename ControlVectorType>
        struct SegmentCalcFirstOrderVisitor
            : fusion::PhaseUnaryVisitorBase<SegmentCalcFirstOrderVisitor<StateVectorType, ControlVectorType>>
        {
            using ArgsType = boost::fusion::vector<std::size_t, Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<ControlVectorType>>;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const std::size_t &segment_index,
                const Eigen::MatrixBase<StateVectorType> &xs,
                const Eigen::MatrixBase<ControlVectorType> &us)
            {
                phase_model.segments_[segment_index].calcDiff(phase_data[segment_index], xs.derived(), us.derived());
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename ControlVectorType>
        inline void segment_calc_first_order(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &segment_index,
            const Eigen::MatrixBase<StateVectorType> &xs,
            const Eigen::MatrixBase<ControlVectorType> &us)
        {
            typedef SegmentCalcFirstOrderVisitor Algo;

            Algo::run(phase_model, phase_data, typename Algo::ArgsType(segment_index, xs, us));
        }

        template <
            typename NumScalar,
            typename StateVectorType,
            typename ControlVectorType>
        struct SegmentQuasiStaticVisitor
            : fusion::PhaseUnaryVisitorBase<SegmentQuasiStaticVisitor<NumScalar, StateVectorType, ControlVectorType>>
        {
            using ArgsType = boost::fusion::vector<std::size_t, Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<ControlVectorType>, const std::size_t &, const NumScalar &>;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const std::size_t &segment_index,
                const Eigen::MatrixBase<StateVectorType> &xs,
                Eigen::MatrixBase<ControlVectorType> &us,
                const std::size_t &maxiter,
                const NumScalar &tol)
            {
                phase_model.segments_[segment_index].quasiStatic(phase_data[segment_index], xs.derived(), us.derived(), maxiter, tol);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename ControlVectorType>
        inline void segment_quasi_static(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            PhaseDataTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_data,
            const std::size_t &segment_index,
            const Eigen::MatrixBase<StateVectorType> &xs,
            Eigen::MatrixBase<ControlVectorType> &us,
            const std::size_t &maxiter,
            const NumScalar &tol)
        {
            typedef SegmentQuasiStaticVisitor<NumScalar, StateVectorType, ControlVectorType> Algo;

            Algo::run(phase_model, phase_data, typename Algo::ArgsType(segment_index, xs, us, maxiter, tol));
        }

        struct PhaseNxVisitor : boost::static_visitor<int>
        {
            template <typename PhaseModelDerived>
            int operator()(const PhaseModelBase<PhaseModelDerived> &phase_model) const
            {
                return phase_model.nx();
            }

            template <typename VarScalar,
                      typename NumScalar,
                      int Options,
                      template <typename, typename, int> class PhaseCollectionTpl>
            static int run(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
            {
                return boost::apply_visitor(PhaseNxVisitor(), phase_model);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int nx(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            return PhaseNxVisitor::run(phase_model);
        }

        struct PhaseNuVisitor : boost::static_visitor<int>
        {
            template <typename PhaseModelDerived>
            int operator()(const PhaseModelBase<PhaseModelDerived> &phase_model) const
            {
                return phase_model.nu();
            }

            template <typename VarScalar,
                      typename NumScalar,
                      int Options,
                      template <typename, typename, int> class PhaseCollectionTpl>
            static int run(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
            {
                return boost::apply_visitor(PhaseNuVisitor(), phase_model);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int nu(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            return PhaseNuVisitor::run(phase_model);
        }

        struct PhaseNdxVisitor : boost::static_visitor<int>
        {
            template <typename PhaseModelDerived>
            int operator()(const PhaseModelBase<PhaseModelDerived> &phase_model) const
            {
                return phase_model.ndx();
            }

            template <typename VarScalar,
                      typename NumScalar,
                      int Options,
                      template <typename, typename, int> class PhaseCollectionTpl>
            static int run(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
            {
                return boost::apply_visitor(PhaseNdxVisitor(), phase_model);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int ndx(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            return PhaseNdxVisitor::run(phase_model);
        }

        struct PhaseNhVisitor : boost::static_visitor<int>
        {
            template <typename PhaseModelDerived>
            int operator()(const PhaseModelBase<PhaseModelDerived> &phase_model) const
            {
                return phase_model.nh();
            }

            template <typename VarScalar,
                      typename NumScalar,
                      int Options,
                      template <typename, typename, int> class PhaseCollectionTpl>
            static int run(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
            {
                return boost::apply_visitor(PhaseNhVisitor(), phase_model);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int nh(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            return PhaseNhVisitor::run(phase_model);
        }

        struct PhaseNgVisitor : boost::static_visitor<int>
        {
            template <typename PhaseModelDerived>
            int operator()(const PhaseModelBase<PhaseModelDerived> &phase_model) const
            {
                return phase_model.ng();
            }

            template <typename VarScalar,
                      typename NumScalar,
                      int Options,
                      template <typename, typename, int> class PhaseCollectionTpl>
            static int run(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
            {
                return boost::apply_visitor(PhaseNgVisitor(), phase_model);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int ng(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            return PhaseNgVisitor::run(phase_model);
        }

        struct PhaseNcVisitor : boost::static_visitor<int>
        {
            template <typename PhaseModelDerived>
            int operator()(const PhaseModelBase<PhaseModelDerived> &phase_model) const
            {
                return phase_model.nc();
            }

            template <typename VarScalar,
                      typename NumScalar,
                      int Options,
                      template <typename, typename, int> class PhaseCollectionTpl>
            static int run(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
            {
                return boost::apply_visitor(PhaseNcVisitor(), phase_model);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline int nc(const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            return PhaseNcVisitor::run(phase_model);
        }

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_basic_visitors_hxx__