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

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        struct PhasePeriodVisitor
            : fusion::PhaseUnaryVisitorBase<PhasePeriodVisitor<VarScalar, NumScalar, Options, PhaseCollectionTpl>>
        {
            using ArgsType = boost::fusion::vector<>;

            template <typename PhaseModel>
            static auto algo(
                const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
            {
                return phase_model.phase_period_;
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline NumScalar phase_period(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            typedef PhasePeriodVisitor<VarScalar, NumScalar, Options, PhaseCollectionTpl> Algo;

            return Algo::run(phase_model, typename Algo::ArgsType());
        }

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::VectorNX_t state_zero(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            typedef StateZeroVisitor<NumScalar> Algo;

            return Algo::run(phase_model, typename Algo::ArgsType());
        }

        template <typename NumScalar>
        struct StateRandVisitor
            : fusion::PhaseUnaryVisitorBase<StateRandVisitor<NumScalar>>
        {
            using ArgsType = boost::fusion::vector<>;

            template <typename PhaseModel>
            static auto algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data)
            {
                return phase_model.state_.rand();
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl>
        inline typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::VectorNX_t state_rand(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model)
        {
            typedef StateRandVisitor<NumScalar> Algo;

            return Algo::run(phase_model, typename Algo::ArgsType());
        }

        template <typename StateVectorType1, typename StateVectorType2, typename StateTangentVectorType>
        struct StateDiffVisitor
            : fusion::PhaseUnaryVisitorBase<StateDiffVisitor<StateVectorType1, StateVectorType2, StateTangentVectorType>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType1>, Eigen::MatrixBase<StateVectorType2>, Eigen::MatrixBase<StateTangentVectorType>>;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType1> &xs,
                const Eigen::MatrixBase<StateVectorType2> &xs_next,
                Eigen::MatrixBase<StateTangentVectorType> &dxout)
            {
                phase_model.state_.diff(xs.derived(), xs_next.derived(), dxout.derived());
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType1,
                  typename StateVectorType2,
                  typename StateTangentVectorType>
        inline void state_diff(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType1> &xs,
            const Eigen::MatrixBase<StateVectorType2> &xs_next,
            Eigen::MatrixBase<StateTangentVectorType> &dxout)
        {
            typedef StateDiffVisitor<StateVectorType1, StateVectorType2, StateTangentVectorType> Algo;

            Algo::run(phase_model, typename Algo::ArgsType(xs, xs_next, dxout));
        }

        template <typename StateVectorType, typename StateTangentVectorType, typename StateVectorType2>
        struct StateIntegrateVisitor
            : fusion::PhaseUnaryVisitorBase<StateIntegrateVisitor<StateVectorType, StateTangentVectorType, StateVectorType2>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<StateTangentVectorType>, Eigen::MatrixBase<StateVectorType2>>;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType> &x,
                const Eigen::MatrixBase<StateTangentVectorType> &dx,
                Eigen::MatrixBase<StateVectorType2> &xout)
            {
                phase_model.state_.integrate(x.derived(), dx.derived(), xout.derived());
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType,
                  typename StateVectorType2>
        inline void state_integrate(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx,
            Eigen::MatrixBase<StateVectorType2> &xout)
        {
            typedef StateIntegrateVisitor<StateVectorType, StateTangentVectorType, StateVectorType2> Algo;

            Algo::run(phase_model, typename Algo::ArgsType(x, dx, xout));
        }

        template <typename StateVectorType1, typename StateVectorType2, typename JMatrix1, typename JMatrix2>
        struct StateJdiffVisitor
            : fusion::PhaseUnaryVisitorBase<StateJdiffVisitor<StateVectorType1, StateVectorType2, JMatrix1, JMatrix2>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType1>, Eigen::MatrixBase<StateVectorType2>, Eigen::MatrixBase<JMatrix1>, Eigen::MatrixBase<JMatrix2>, galileo::core::Jcomponent>;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType1> &x0,
                const Eigen::MatrixBase<StateVectorType2> &x1,
                Eigen::MatrixBase<JMatrix1> &Jfirst,
                Eigen::MatrixBase<JMatrix2> &Jsecond,
                const galileo::core::Jcomponent firstsecond)
            {
                phase_model.state_.Jdiff(x0.derived(), x1.derived(), Jfirst.derived(), Jsecond.derived(), firstsecond);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType1,
                  typename StateVectorType2,
                  typename JMatrix1,
                  typename JMatrix2>
        inline void state_jdiff(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType1> &x0,
            const Eigen::MatrixBase<StateVectorType2> &x1,
            Eigen::MatrixBase<JMatrix1> &Jfirst,
            Eigen::MatrixBase<JMatrix2> &Jsecond,
            const galileo::core::Jcomponent firstsecond)
        {
            typedef StateJdiffVisitor<StateVectorType1, StateVectorType2, JMatrix1, JMatrix2> Algo;

            Algo::run(phase_model, typename Algo::ArgsType(x0, x1, Jfirst, Jsecond, firstsecond));
        }

        template <typename StateVectorType, typename StateTangentVectorType, typename JMatrix1, typename JMatrix2>
        struct StateJintegrateVisitor
            : fusion::PhaseUnaryVisitorBase<StateJintegrateVisitor<StateVectorType, StateTangentVectorType, JMatrix1, JMatrix2>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<StateTangentVectorType>, Eigen::MatrixBase<JMatrix1>, Eigen::MatrixBase<JMatrix2>, galileo::core::Jcomponent, galileo::core::AssignmentOp>;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType> &x,
                const Eigen::MatrixBase<StateTangentVectorType> &dx,
                Eigen::MatrixBase<JMatrix1> &Jfirst,
                Eigen::MatrixBase<JMatrix2> &Jsecond,
                const galileo::core::Jcomponent firstsecond,
                const galileo::core::AssignmentOp op)
            {
                phase_model.state_.Jintegrate(x.derived(), dx.derived(), Jfirst.derived(), Jsecond.derived(), firstsecond, op);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType,
                  typename JMatrix1,
                  typename JMatrix2>
        inline void state_jintegrate(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx,
            Eigen::MatrixBase<JMatrix1> &Jfirst,
            Eigen::MatrixBase<JMatrix2> &Jsecond,
            const galileo::core::Jcomponent firstsecond,
            const galileo::core::AssignmentOp op)
        {
            typedef StateJintegrateVisitor<StateVectorType, StateTangentVectorType, JMatrix1, JMatrix2> Algo;

            Algo::run(phase_model, typename Algo::ArgsType(x, dx, Jfirst, Jsecond, firstsecond, op));
        }

        template <typename StateVectorType, typename StateTangentVectorType, typename JMatrix>
        struct StateJintegrateTransportVisitor
            : fusion::PhaseUnaryVisitorBase<StateJintegrateTransportVisitor<StateVectorType, StateTangentVectorType, JMatrix>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<StateTangentVectorType>, Eigen::MatrixBase<JMatrix>, galileo::core::Jcomponent>;

            template <typename PhaseModel>
            static void algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType> &x,
                const Eigen::MatrixBase<StateTangentVectorType> &dx,
                Eigen::MatrixBase<JMatrix> &Jin,
                const galileo::core::Jcomponent firstsecond)
            {
                phase_model.state_.JintegrateTransport(x.derived(), dx.derived(), Jin.derived(), firstsecond);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType,
                  typename JMatrix>
        inline void state_jintegrate_transport(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx,
            Eigen::MatrixBase<JMatrix> &Jin,
            const galileo::core::Jcomponent firstsecond)
        {
            typedef StateJintegrateTransportVisitor<StateVectorType, StateTangentVectorType, JMatrix> Algo;

            Algo::run(phase_model, typename Algo::ArgsType(x, dx, Jin, firstsecond));
        }

        template <typename StateVectorType1, typename StateVectorType2>
        struct StateDiffDxVisitor
            : fusion::PhaseUnaryVisitorBase<StateDiffDxVisitor<StateVectorType1, StateVectorType2>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType1>, Eigen::MatrixBase<StateVectorType2>>;

            template <typename PhaseModel>
            static auto algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType1> &x0,
                const Eigen::MatrixBase<StateVectorType2> &x1)
            {
                return phase_model.state_.diff_dx(x0.derived(), x1.derived());
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType>
        inline typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::VectorNX_t state_diff_dx(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x0,
            const Eigen::MatrixBase<StateVectorType> &x1)
        {
            typedef StateDiffDxVisitor<StateVectorType, StateVectorType> Algo;

            return Algo::run(phase_model, typename Algo::ArgsType(x0, x1));
        }

        template <typename StateVectorType, typename StateTangentVectorType>
        struct StateIntegrateXVisitor
            : fusion::PhaseUnaryVisitorBase<StateIntegrateXVisitor<StateVectorType, StateTangentVectorType>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<StateTangentVectorType>>;

            template <typename PhaseModel>
            static auto algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType> &x,
                const Eigen::MatrixBase<StateTangentVectorType> &dx)
            {
                return phase_model.state_.integrate_x(x.derived(), dx.derived());
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType>
        inline typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::VectorNX_t state_integrate_x(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx)
        {
            typedef StateIntegrateXVisitor<StateVectorType, StateTangentVectorType> Algo;

            return Algo::run(phase_model, typename Algo::ArgsType(x, dx));
        }

        template <typename StateVectorType1, typename StateVectorType2>
        struct StateJdiffJsVisitor
            : fusion::PhaseUnaryVisitorBase<StateJdiffJsVisitor<StateVectorType1, StateVectorType2>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType1>, Eigen::MatrixBase<StateVectorType2>, galileo::core::Jcomponent>;

            template <typename PhaseModel>
            static auto algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType1> &x0,
                const Eigen::MatrixBase<StateVectorType2> &x1,
                const galileo::core::Jcomponent firstsecond)
            {
                return phase_model.state_.Jdiff_Js(x0.derived(), x1.derived(), firstsecond);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType1,
                  typename StateVectorType2>
        inline std::vector<typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::MatrixNDX_t> state_jdiff_Js(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType1> &x0,
            const Eigen::MatrixBase<StateVectorType2> &x1,
            const galileo::core::Jcomponent firstsecond)
        {
            typedef StateJdiffJsVisitor<StateVectorType1, StateVectorType2> Algo;

            return Algo::run(phase_model, typename Algo::ArgsType(x0, x1, firstsecond));
        }

        template <typename StateVectorType, typename StateTangentVectorType>
        struct StateJintegrateJsVisitor
            : fusion::PhaseUnaryVisitorBase<StateJintegrateJsVisitor<StateVectorType, StateTangentVectorType>>
        {
            using ArgsType = boost::fusion::vector<Eigen::MatrixBase<StateVectorType>, Eigen::MatrixBase<StateTangentVectorType>, galileo::core::Jcomponent>;

            template <typename PhaseModel>
            static auto algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data,
                const Eigen::MatrixBase<StateVectorType> &x,
                const Eigen::MatrixBase<StateTangentVectorType> &dx,
                const galileo::core::Jcomponent firstsecond)
            {
                return phase_model.state_.Jintegrate_Js(x.derived(), dx.derived(), firstsecond);
            }
        };

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  template <typename, typename, int> class PhaseCollectionTpl,
                  typename StateVectorType,
                  typename StateTangentVectorType>
        inline std::vector<typename PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl>::State_t::MatrixNDX_t> state_jintegrate_Js(
            const PhaseModelTpl<VarScalar, NumScalar, Options, PhaseCollectionTpl> &phase_model,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<StateTangentVectorType> &dx,
            const galileo::core::Jcomponent firstsecond)
        {
            typedef StateJintegrateJsVisitor<StateVectorType, StateTangentVectorType> Algo;

            return Algo::run(phase_model, typename Algo::ArgsType(x, dx, firstsecond));
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

        template <typename NumScalar>
        struct StateZeroVisitor
            : fusion::PhaseUnaryVisitorBase<StateZeroVisitor<NumScalar>>
        {
            using ArgsType = boost::fusion::vector<>;

            template <typename PhaseModel>
            static auto algo(
                const galileo::predictive::PhaseModelBase<PhaseModel> &phase_model,
                typename galileo::predictive::PhaseDataBase<typename PhaseModel::PhaseDataDerived> &phase_data)
            {
                return phase_model.state_.zero();
            }
        };

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_basic_visitors_hxx__