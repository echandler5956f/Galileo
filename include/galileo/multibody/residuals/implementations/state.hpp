#ifndef __galileo_multibody_residuals_state_hpp__
#define __galileo_multibody_residuals_state_hpp__

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/core/states/state-base.hpp"
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

        using DimNR_t = DimensionTpl<PS::DimNDX_t::Value>;
        static constexpr int NR = DimNR_t::Value;

        static constexpr bool QDependent = true;
        static constexpr bool VDependent = true;
        static constexpr bool UDependent = false;

        using R_t = Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>;
        using Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::DimNDX_t::Value, PS::Options>;
        using Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::DimNU_t::Value, PS::Options>;
        using Arr_Rx_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::DimNDX_t::Value, PS::Options>;
        using Arr_Ru_t = Eigen::GMatrix<typename PS::VarScalar, NR, PS::DimNU_t::Value, PS::Options>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualDataStateTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualStateTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelStateTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualStateTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct ResidualDataStateTpl
        : public ResidualDataBase<ResidualDataStateTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

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

        ResidualDataStateTpl(const Model_t &model)
            : R(model.get_nr()), Rx(model.get_nr(), model.get_ps().ndx_dim.value()),
              Ru(model.get_nr(), model.get_ps().nu_dim.value()),
              Arr_Rx(model.get_nr(), model.get_ps().ndx_dim.value()),
              Arr_Ru(model.get_nr(), model.get_ps().nu_dim.value())
        {
            R.setZero();
            Rx.setZero();
            Ru.setZero();
            Arr_Rx.setZero();
            Arr_Ru.setZero();
        }

        typename PS::RobotData_t *robot;

        R_t R;
        Rx_t Rx;
        Ru_t Ru;
        Arr_Rx_t Arr_Rx;
        Arr_Ru_t Arr_Ru;

    }; // class ResidualDataStateTpl

    template <typename PhaseSpec>
    class ResidualModelStateTpl
        : public ResidualModelBase<ResidualModelStateTpl<PhaseSpec>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = ResidualStateTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelBase<ResidualModelStateTpl<PS>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        ResidualModelStateTpl(const PS &ps,
                              const std::shared_ptr<State_t> &state,
                              const VectorNx_t &x_ref)
            : Base(ps, DimNR_t()),
              state_(state), x_ref_(x_ref)
        {
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            state_->diff(x_ref_, x, data.R);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            state_->Jdiff(x_ref_, x, data.Rx, data.Rx, Jcomponent::second);
        }

        template <typename CostDataType, typename ActivationDataType, bool UpdateU = true>
        void calcCostDiffImpl(CostDataType &cdata,
                              Data_t &rdata,
                              const ActivationDataType &adata) const
        {
            const PS &ps = get_ps();
            const RobotModel_t &robot = state_->get_robot();
            typedef Eigen::Block<MatrixX_t> MatrixBlock;

            // trust
            for (pinocchio::JointIndex i = 1;
                 i < (pinocchio::JointIndex)robot.njoints; ++i)
            {
                const MatrixBlock &RxBlock = block(rdata.Rx, robot.idx_vs[i], robot.idx_vs[i],
                                                   robot.nvs[i], robot.nvs[i]);
                segment(cdata.Lx, robot.idx_vs[i], robot.nvs[i]).noalias() =
                    RxBlock.transpose() *
                    segment(adata.Ar, robot.idx_vs[i], robot.nvs[i]);

                block(cdata.Lxx, robot.idx_vs[i], robot.idx_vs[i], robot.nvs[i], robot.nvs[i])
                    .noalias() = RxBlock.transpose() *
                                 segment(adata.Arr.diagonal(), robot.idx_vs[i], robot.nvs[i])
                                     .asDiagonal() *
                                 RxBlock;
            }
            tail(cdata.Lx, ps.nv_dim) = tail(adata.Ar, ps.nv_dim);
            tail(cdata.Lxx.diagonal(), ps.nv_dim) =
                tail(adata.Arr.diagonal(), ps.nv_dim);
        }

        template <typename DataCollector>
        Data_t createData(DataCollector *const collector) const
        {
            Data_t data(*this);
            data.robot = collector->robot;
            return data;
        }

        const std::shared_ptr<State_t> &get_state() const
        {
            return state_;
        }

        using Base::get_ps;

        using Base::get_nr;
        using Base::get_nr_dim;

        using Base::get_nu;
        using Base::get_nu_dim;

        using Base::get_q_dependent;
        using Base::get_v_dependent;
        using Base::get_u_dependent;

    protected:
        std::shared_ptr<State_t> state_;
        VectorNx_t x_ref_;

    }; // class ResidualModelStateTpl

} // namespace galileo

#endif // __galileo_multibody_residuals_state_hpp__
