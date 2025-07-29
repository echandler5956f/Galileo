#ifndef __galileo_predictive_phases_phase_multibody_hpp__
#define __galileo_predictive_phases_phase_multibody_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <
        typename PhaseSpec,
        template <typename _PS, template <typename PS_> class ResetCollectionTpl> class ResetNodeTpl>
    struct PhaseMultibodyTpl;

    template <typename PhaseSpec,
              template <typename _PS, template <typename PS_> class ResetCollectionTpl> class ResetNodeTpl>
    struct traits<PhaseMultibodyTpl<PhaseSpec, ResetNodeTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseMultibodyTpl<PS, ResetNodeTpl>;
        using Model_t = PhaseModelMultibodyTpl<PS, ResetNodeTpl>;
        using Data_t = PhaseDataMultibodyTpl<PS, ResetNodeTpl>;

        using ResetNodeMeta_t = typename traits<ResetNodeTpl<PS, ResetCollectionTpl>>::Meta_t;
        using ResetNodeModel_t = typename traits<ResetNodeMeta_t>::Model_t;
        using ResetNodeData_t = typename traits<ResetNodeMeta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename _PS, template <typename PS_> class ResetCollectionTpl> class ResetNodeTpl>
    struct traits<PhaseDataMultibodyTpl<PhaseSpec, ResetNodeTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseMultibodyTpl<PS, ResetNodeTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename _PS, template <typename PS_> class ResetCollectionTpl> class ResetNodeTpl>
    struct traits<PhaseModelMultibodyTpl<PhaseSpec, ResetNodeTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseMultibodyTpl<PS, ResetNodeTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename _PS, template <typename PS_> class ResetCollectionTpl> class ResetNodeTpl>
    struct PhaseDataMultibodyTpl
        : public PhaseDataBase<PhaseDataMultibodyTpl<PhaseSpec, ResetNodeTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseMultibodyTpl<PS, ResetNodeTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseDataBase<PhaseDataMultibodyTpl<PS, ResetNodeTpl>, PS>;

        using ResetNodeMeta_t = typename traits<Meta_t>::ResetNodeMeta_t;
        using ResetNodeModel_t = typename traits<ResetNodeMeta_t>::Model_t;
        using ResetNodeData_t = typename traits<ResetNodeMeta_t>::Data_t;

        PhaseDataMultibodyTpl(const Model_t &model)
            : Base(model)
        {
        }

        using Base::get_segments;

        using Base::XNext_at_i;
        using Base::XNextw_at_i;
        using Base::XNextx_at_i;

        using Base::L_at_i;
        using Base::Lw_at_i;
        using Base::Lx_at_i;

        using Base::Lww_at_i;
        using Base::Lxw_at_i;
        using Base::Lxx_at_i;

        using Base::H_at_i;
        using Base::Hw_at_i;
        using Base::Hx_at_i;

        using Base::G_at_i;
        using Base::Gw_at_i;
        using Base::Gx_at_i;

        ResetNodeData_t reset_node_data;

    }; // struct PhaseDataMultibodyTpl

    template <typename PhaseSpec,
              template <typename _PS, template <typename PS_> class ResetCollectionTpl> class ResetNodeTpl>
    class PhaseModelMultibodyTpl
        : public PhaseModelBase<PhaseModelMultibodyTpl<PhaseSpec, ResetNodeTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseMultibodyTpl<PS, ResetNodeTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseModelBase<PhaseModelMultibodyTpl<PS, ResetNodeTpl>, PS>;

        using ResetNodeMeta_t = typename traits<Meta_t>::ResetNodeMeta_t;
        using ResetNodeModel_t = typename traits<ResetNodeMeta_t>::Model_t;
        using ResetNodeData_t = typename traits<ResetNodeMeta_t>::Data_t;

        PhaseModelMultibodyTpl(const PS &ps)
            : Base(ps)
        {
        }

        template <bool IsForward, typename RightPhaseModelType, typename RightPhaseDataType, typename StateVectorType, typename ControlVectorType>
        void resetMapCalc(Data_t &data,
                          const RightPhaseModelType &right_model,
                          RightPhaseDataType &right_data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            if constexpr (IsForward)
            {
                right_model.calc(right_data, x, u);
            }
            else
            {
                reset_node_.calc(data.reset_node_data, x, u);
            }
        }

        template <bool IsForward, typename RightPhaseModelType, typename RightPhaseDataType, typename StateVectorType>
        void resetMapCalc(Data_t &data,
                          const RightPhaseModelType &right_model,
                          RightPhaseDataType &right_data,
                          const Eigen::MatrixBase<StateVectorType> &x) const
        {
            if constexpr (IsForward)
            {
                right_model.calc(right_data, x);
            }
            else
            {
                reset_node_.calc(data.reset_node_data, x);
            }
        }

        template <bool IsForward, typename RightPhaseModelType, typename RightPhaseDataType, typename StateVectorType, typename ControlVectorType>
        void resetMapCalcDiff(Data_t &data,
                              const RightPhaseModelType &right_model,
                              RightPhaseDataType &right_data,
                              const Eigen::MatrixBase<StateVectorType> &x,
                              const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            if constexpr (IsForward)
            {
                right_model.calcDiff(right_data, x, u);
            }
            else
            {
                reset_node_.calcDiff(data.reset_node_data, x, u);
            }
        }

        template <bool IsForward, typename RightPhaseModelType, typename RightPhaseDataType, typename StateVectorType>
        void resetMapCalcDiff(Data_t &data,
                              const RightPhaseModelType &right_model,
                              RightPhaseDataType &right_data,
                              const Eigen::MatrixBase<StateVectorType> &x) const
        {
            if constexpr (IsForward)
            {
                right_model.calcDiff(right_data, x);
            }
            else
            {
                reset_node_.calcDiff(data.reset_node_data, x);
            }
        }

        using Base::calc;
        using Base::calcDiff;
        using Base::quasiStatic;

        using Base::createData;

        using Base::get_ps;
        using Base::get_segments;

    protected:
        ResetNodeModel_t reset_node_;
    }; // class PhaseModelMultibodyTpl

} // namespace galileo

#endif // __galileo_predictive_phases_phase_multibody_hpp__
