#ifndef __galileo_multibody_robot_holder_hpp__
#define __galileo_multibody_robot_holder_hpp__

#include "galileo/core/fwd.hpp"

namespace galileo
{

    template <int _NQb,
              int _NQj,
              int _NVb,
              int _NVj,
              int _NRotors>
    struct RobotHolderTpl
    {
        using DimNQb_t = DimensionTpl<_NQb>;                    // Dimension of floating base generalized coordinates
        using DimNQj_t = DimensionTpl<_NQj>;                    // Dimension of joint generalized coordinates
        using DimNVb_t = DimensionTpl<_NVb>;                    // Dimension of floating base generalized velocities
        using DimNVj_t = DimensionTpl<_NVj>;                    // Dimension of joint generalized velocities
        using DimNRotors_t = DimensionTpl<_NRotors>;            // Number of rotors attached to the floating base
        using DimNQ_t = decltype(DimNQb_t{} + DimNQj_t{});      // Dimension of generalized coordinates
        using DimNV_t = decltype(DimNVb_t{} + DimNVj_t{});      // Dimension of generalized velocities
        using DimNX_t = decltype(DimNQ_t{} + DimNV_t{});        // Dimension of state
        using DimNDX_t = decltype(DimNV_t{} + DimNV_t{});       // Dimension of state tangent space
        using DimNUa_t = decltype(DimNVj_t{} + DimNRotors_t{}); // Dimension of actuated inputs

        static constexpr int NQb = DimNQb_t::Value;
        static constexpr int NQj = DimNQj_t::Value;
        static constexpr int NVb = DimNVb_t::Value;
        static constexpr int NVj = DimNVj_t::Value;
        static constexpr int NRotors = DimNRotors_t::Value;
        static constexpr int NQ = DimNQ_t::Value;
        static constexpr int NV = DimNV_t::Value;
        static constexpr int NX = DimNX_t::Value;
        static constexpr int NDX = DimNDX_t::Value;
        static constexpr int NUa = DimNUa_t::Value;

        RobotHolderTpl()
            : nqb_dim_{}, nqj_dim_{}, nvb_dim_{}, nvj_dim_{}, nrotors_dim_{},
              nq_dim_(nqb_dim_ + nqj_dim_),
              nv_dim_(nvb_dim_ + nvj_dim_),
              nx_dim_(nq_dim_ + nv_dim_),
              ndx_dim_(nv_dim_ + nv_dim_),
              nua_dim_(nvj_dim_ + nrotors_dim_)
        {
        }

        const DimNQb_t &get_nqb_dim() const
        {
            return nqb_dim_;
        }

        int get_nqb() const
        {
            return nqb_dim_.value();
        }

        const DimNQj_t &get_nqj_dim() const
        {
            return nqj_dim_;
        }

        int get_nqj() const
        {
            return nqj_dim_.value();
        }

        const DimNVb_t &get_nvb_dim() const
        {
            return nvb_dim_;
        }

        int get_nvb() const
        {
            return nvb_dim_.value();
        }

        const DimNVj_t &get_nvj_dim() const
        {
            return nvj_dim_;
        }

        int get_nvj() const
        {
            return nvj_dim_.value();
        }

        const DimNRotors_t &get_nrotors_dim() const
        {
            return nrotors_dim_;
        }

        int get_nrotors() const
        {
            return nrotors_dim_.value();
        }

        const DimNQ_t &get_nq_dim() const
        {
            return nq_dim_;
        }

        int get_nq() const
        {
            return nq_dim_.value();
        }

        const DimNV_t &get_nv_dim() const
        {
            return nv_dim_;
        }

        int get_nv() const
        {
            return nv_dim_.value();
        }

        const DimNX_t &get_nx_dim() const
        {
            return nx_dim_;
        }

        int get_nx() const
        {
            return nx_dim_.value();
        }

        const DimNDX_t &get_ndx_dim() const
        {
            return ndx_dim_;
        }

        int get_ndx() const
        {
            return ndx_dim_.value();
        }

        const DimNUa_t &get_nua_dim() const
        {
            return nua_dim_;
        }

        int get_nua() const
        {
            return nua_dim_.value();
        }

        DimNQb_t nqb_dim_;
        DimNQj_t nqj_dim_;
        DimNVb_t nvb_dim_;
        DimNVj_t nvj_dim_;
        DimNRotors_t nrotors_dim_;
        DimNQ_t nq_dim_;
        DimNV_t nv_dim_;
        DimNX_t nx_dim_;
        DimNDX_t ndx_dim_;
        DimNUa_t nua_dim_;

        void display(std::ostream &os, const std::string &indent = "  ") const
        {
            os << indent << "Robot Dimensions:\n";
            os << indent << "  Base Configuration:\n";
            os << indent << "    NQb (Base generalized coordinates): " << get_nqb_dim() << "\n";
            os << indent << "    NVb (Base generalized velocities):  " << get_nvb_dim() << "\n";
            os << indent << "  Joint Configuration:\n";
            os << indent << "    NQj (Joint generalized coordinates): " << get_nqj_dim() << "\n";
            os << indent << "    NVj (Joint generalized velocities):  " << get_nvj_dim() << "\n";
            os << indent << "  Actuation:\n";
            os << indent << "    NRotors (Number of rotors): " << get_nrotors_dim() << "\n";
            os << indent << "    NUa (Actuated inputs):      " << get_nua_dim() << "\n";
            os << indent << "  Combined Dimensions:\n";
            os << indent << "    NQ (Total generalized coordinates): " << get_nq_dim() << "\n";
            os << indent << "    NV (Total generalized velocities):  " << get_nv_dim() << "\n";
            os << indent << "    NX (State dimension):                " << get_nx_dim() << "\n";
            os << indent << "    NDX (State tangent dimension):       " << get_ndx_dim() << "\n";
            os << indent << "Configuration: " << (IsValidRobotHolder(*this) ? "VALID" : "INVALID") << "\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const RobotHolderTpl &rh)
        {
            os << "RobotHolder{\n";
            rh.display(os, "  ");
            os << "}";
            return os;
        }
    };

    template <typename RobotHolderType>
    bool IsValidRobotHolder(const RobotHolderType &rh)
    {
        bool valid_nq = (rh.get_nq() == rh.get_nqb() + rh.get_nqj());
        bool valid_nv = (rh.get_nv() == rh.get_nvb() + rh.get_nvj());
        bool valid_nx = (rh.get_nx() == rh.get_nq() + rh.get_nv());
        bool valid_ndx = (rh.get_ndx() == rh.get_nv() + rh.get_nv());
        bool valid_nua = (rh.get_nua() == rh.get_nvj() + rh.get_nrotors());
        return valid_nq && valid_nv && valid_nx && valid_ndx && valid_nua;
    }

} // namespace galileo

#endif // __galileo_multibody_robot_holder_hpp__
