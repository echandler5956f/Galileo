#ifndef __galileo_predictive_phases_phase_data_base_hpp__
#define __galileo_predictive_phases_phase_data_base_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-model-base.hpp"

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct PhaseDataBase : internal::CRTP<PhaseDataBase<Derived, PhaseSpec>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using SegmentDataVector_t = typename PS::SegmentDataVector_t;

        FORWARD_ACCESSOR(SegmentDataVector_t, segments);

        // The fully expanded contents of each PhaseData derived class should be
        // SegmentDataVector segments;
        // ---For each segment data, we have---
        // | - C_t C;
        // | - Ck_t Ck;
        // | - Cw_t Cw;
        // | - ControlParamData_t control;
        // |   | - W_t W;
        // |   | - U_t U;
        // |   | - Wu_t Wu;
        // | - ConstraintDataCollection_t constraints;
        // |   | - ResidualData_t residual;
        // |   |   | - R_t R;
        // |   |   | - Rx_t Rx;
        // |   |   | - Ru_t Ru;
        // |   |   | - Arr_Rx_t Arr_Rx;
        // |   |   | - Arr_Ru_t Arr_Ru;
        // |   | - Seg_H_t H;
        // |   | - Seg_Hx_t Hx;
        // |   | - Seg_Hu_t Hu;
        // |   | - Seg_G_t G;
        // |   | - Seg_Gx_t Gx;
        // |   | - Seg_Gu_t Gu;
        // | - CostDataCollection_t costs;
        // |   | - ResidualData_t residual;
        // |   |   | - R_t R;
        // |   |   | - Rx_t Rx;
        // |   |   | - Ru_t Ru;
        // |   |   | - Arr_Rx_t Arr_Rx;
        // |   |   | - Arr_Ru_t Arr_Ru;
        // |   | - ActivationData_t activation;
        // |   |   | - A_t A;
        // |   |   | - Ar_t Ar;
        // |   |   | - Arr_t Arr;
        // |   | - L_t L;
        // |   | - Lx_t Lx;
        // |   | - Lu_t Lu;
        // |   | - Lxx_t Lxx;
        // |   | - Lxu_t Lxu;
        // |   | - Luu_t Luu;
        // | - NodeDataVector nodes;
        // | ---For each node data, we have---
        // |   | - ActuationData_t actuation;
        // |   | - ConstraintDataCollection_t constraints;
        // |   |   | - ResidualData_t residual;
        // |   |   |   | - R_t R;
        // |   |   |   | - Rx_t Rx;
        // |   |   |   | - Ru_t Ru;
        // |   |   |   | - Arr_Rx_t Arr_Rx;
        // |   |   |   | - Arr_Ru_t Arr_Ru;
        // |   |   | - H_t H;
        // |   |   | - Hx_t Hx;
        // |   |   | - Hu_t Hu;
        // |   |   | - G_t G;
        // |   |   | - Gx_t Gx;
        // |   |   | - Gu_t Gu;
        // |   | - CostDataCollection_t costs;
        // |   |   | - ResidualData_t residual;
        // |   |   |   | - R_t R;
        // |   |   |   | - Rx_t Rx;
        // |   |   |   | - Ru_t Ru;
        // |   |   |   | - Arr_Rx_t Arr_Rx;
        // |   |   |   | - Arr_Ru_t Arr_Ru;
        // |   |   | - ActivationData_t activation;
        // |   |   |   | - A_t A;
        // |   |   |   | - Ar_t Ar;
        // |   |   |   | - Arr_t Arr;
        // |   |   | - L_t L;
        // |   |   | - Lx_t Lx;
        // |   |   | - Lu_t Lu;
        // |   |   | - Lxx_t Lxx;
        // |   |   | - Lxu_t Lxu;
        // |   |   | - Luu_t Luu;
        // |   | - DynamicsData_t dynamics;
        // |   |   | - F_t F;
        // |   |   | - Fx_t Fx;
        // |   |   | - Fu_t Fu;

        // Thus, we need generic accessors for each of the above.

        // C_t segment_C(const std::size_t &segment_index) const
        // {
        //     return derived().segment_C(segment_index);
        // }

        // Ck_t segment_Ck(const std::size_t &segment_index) const
        // {
        //     return derived().segment_Ck(segment_index);
        // }

        // Cw_t segment_Cw(const std::size_t &segment_index) const
        // {
        //     return derived().segment_Cw(segment_index);
        // }

    protected:
        inline PhaseDataBase()
        {
        }

        inline PhaseDataBase(const PhaseDataBase &clone)
        {
            *this = clone;
        }

        inline PhaseDataBase &operator=(const PhaseDataBase &clone)
        {
            return *this;
        }

    }; // struct PhaseDataBase

} // namespace galileo

#endif // __galileo_predictive_phases_phase_data_base_hpp__
