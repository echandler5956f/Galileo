#ifndef __galileo_multibody_core_residuals_residual_multibody_state_hpp__
#define __galileo_multibody_core_residuals_residual_multibody_state_hpp__

#include <pinocchio/spatial/motion.hpp>

#include <pinocchio/algorithm/joint-configuration.hpp>

#include "galileo/core/residuals/residual-base.hpp"
#include "galileo/core/residuals/impl/residual-state.hpp"
#include "galileo/domains/multibody/core/residuals/fwd.hpp"

namespace galileo
{
    template <typename PhaseSpec>
    struct ResidualMultibodyStateTpl;

    template <typename PhaseSpec>
    struct traits<ResidualMultibodyStateTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = ResidualMultibodyStateTpl<PS>;
        using Model_t = ResidualModelMultibodyStateTpl<PS>;
        using Data_t = ResidualDataStateTpl<PS>;

        using DimNR_t = typename PS::DimNDX_t;
        static constexpr int NR = DimNR_t::Value;

        static constexpr bool QDependent = true;
        static constexpr bool VDependent = true;
        static constexpr bool UDependent = false;

        using R_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, 1, PS::Options>>;
        using Rx_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>>;
        using Ru_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>>;
        using Arr_Rx_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NDX, PS::Options>>;
        using Arr_Ru_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NR, PS::NU, PS::Options>>;
    };

    template <typename PhaseSpec>
    struct traits<ResidualModelMultibodyStateTpl<PhaseSpec>>
    {
        using Meta_t = ResidualMultibodyStateTpl<PhaseSpec>;
    };

    template <typename PhaseSpec>
    class ResidualModelMultibodyStateTpl : public ResidualModelStateTpl<PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ResidualMultibodyStateTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ResidualModelStateTpl<PS>;

        using State_t = typename PS::State_t;
        using RobotModel_t = typename PS::RobotModel_t;
        using JointIndex_t = typename PS::JointIndex_t;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        template <typename StateVectorType>
        ResidualModelMultibodyStateTpl(const PS &ps,
                                       const State_t &state,
                                       const Eigen::MatrixBase<StateVectorType> &x_ref)
            : Base(ps, state, DimNR_t(ps.get_ndx()), x_ref)
        {
        }

        using Base::calc;
        using Base::calcDiff;

        template <bool UpdateU = true, typename CostDataType, typename ActivationDataType>
        void calcCostDiffImpl(CostDataType &cdata, Data_t &rdata, const ActivationDataType &adata) const
        {
            const PS &ps = get_ps();
            const RobotModel_t &robot = get_state().get_robot();

            for (JointIndex_t i = 1; i < (JointIndex_t) robot.njoints; ++i)
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

        using Base::createData;

        using Base::get_ps;
        using Base::get_state;
        using Base::get_nr;
        using Base::get_nr_dim;
        using Base::get_q_dependent;
        using Base::get_u_dependent;
        using Base::get_v_dependent;

    }; // class ResidualModelMultibodyStateTpl

} // namespace galileo

#endif // __galileo_multibody_core_residuals_residual_multibody_state_hpp__
