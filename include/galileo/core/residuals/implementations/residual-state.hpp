#ifndef __galileo_core_residuals_residual_state_hpp__
#define __galileo_core_residuals_residual_state_hpp__

#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include "galileo/core/states/state-base.hpp"

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct ResidualStateTpl;

    template <typename PhaseSpec>
    struct traits<ResidualStateTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualStateTpl<PS>;
        using Model_t = ResidualModelStateTpl<PS>;
        using Data_t = ResidualDataStateTpl<PS>;

        using DimNR_t = typename PS::DimNDX_t;

        static constexpr bool QDependent = true;
        static constexpr bool VDependent = true;
        static constexpr bool UDependent = false;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNU_t::Value, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNDX_t::Value, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, PS::DimNU_t::Value, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataStateTpl<PhaseSpec>>
    {
        using Meta_t = ResidualStateTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelStateTpl<PhaseSpec>>
    {
        using Meta_t = ResidualStateTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    struct ResidualDataStateTpl : public ResidualDataBase<ResidualDataStateTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = ResidualStateTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualDataBase<ResidualDataStateTpl<PS>, PS>;

        GALILEO_RESIDUAL_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(R_t, R);
        DEFAULT_ACCESSOR(Rx_t, Rx);
        DEFAULT_ACCESSOR(Ru_t, Ru);
        DEFAULT_ACCESSOR(Arr_Rx_t, Arr_Rx);
        DEFAULT_ACCESSOR(Arr_Ru_t, Arr_Ru);

        template <typename DataCollector>
        ResidualDataStateTpl(const Model_t &model, DataCollector *const collector)
            : robot(collector->robot),
              R(model.get_nr()),
              Rx(model.get_nr(), model.get_ps().get_ndx()),
              Ru(model.get_nr(), model.get_ps().get_nu()),
              Arr_Rx(model.get_nr(), model.get_ps().get_ndx()),
              Arr_Ru(model.get_nr(), model.get_ps().get_nu())
        {
            R.setZero();
            Rx.setZero();
            Ru.setZero();
            Arr_Rx.setZero();
            Arr_Ru.setZero();
        }

        RobotData_t *robot;

        R_t R;
        Rx_t Rx;
        Ru_t Ru;
        Arr_Rx_t Arr_Rx;
        Arr_Ru_t Arr_Ru;

    }; // class ResidualDataStateTpl

    template <typename PhaseSpec>
    class ResidualModelStateTpl : public ResidualModelBase<ResidualModelStateTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = ResidualStateTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelStateTpl<PS>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        template <typename StateVectorType>
        ResidualModelStateTpl(const PS &ps, const Eigen::MatrixBase<StateVectorType> &x_ref)
            : Base(ps, DimNR_t()), x_ref_(x_ref)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            get_ps().get_state().diff(x_ref_, x, data.R);
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            calc(data, x, VectorNu_t::Zero(get_ps().get_nu()));
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            get_ps().get_state().template Jdiff<SECOND>(x_ref_, x, data.Rx, data.Rx);
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            calcDiff(data, x, VectorNu_t::Zero(get_ps().get_nu()));
        }

        template <bool UpdateU = true, typename CostDataType, typename ActivationDataType>
        void calcCostDiffImpl(CostDataType &cdata, Data_t &rdata, const ActivationDataType &adata) const
        {
            const PS &ps = get_ps();
            const RobotModel_t &robot = get_ps().get_state().get_robot();

            // trust
            for (pinocchio::JointIndex i = 1; i < (pinocchio::JointIndex) robot.njoints; ++i)
            {
                const auto &RxBlock = block(rdata.Rx, robot.idx_vs[i], robot.idx_vs[i], robot.nvs[i], robot.nvs[i]);
                segment(cdata.Lx, robot.idx_vs[i], robot.nvs[i]).noalias() =
                    RxBlock.transpose() * segment(adata.Ar, robot.idx_vs[i], robot.nvs[i]);

                block(cdata.Lxx, robot.idx_vs[i], robot.idx_vs[i], robot.nvs[i], robot.nvs[i]).noalias() =
                    RxBlock.transpose() * segment(adata.Arr.diagonal(), robot.idx_vs[i], robot.nvs[i]).asDiagonal() *
                    RxBlock;
            }
            tail(cdata.Lx, ps.get_nv_dim()) = tail(adata.Ar, ps.get_nv_dim());
            cdata.Lxx.diagonal().tail(ps.get_nv_dim()).noalias() = adata.Arr.diagonal().tail(ps.get_nv_dim());
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            return Data_t(*this, collector);
        }

        using Base::get_ps;
        using Base::get_nr;
        using Base::get_nr_dim;
        using Base::get_q_dependent;
        using Base::get_u_dependent;
        using Base::get_v_dependent;

    protected:
        VectorNx_t x_ref_;

    }; // class ResidualModelStateTpl

} // namespace galileo

#endif // __galileo_core_residuals_residual_state_hpp__
