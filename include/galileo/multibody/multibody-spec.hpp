#ifndef __galileo_multibody_multibody_spec_hpp__
#define __galileo_multibody_multibody_spec_hpp__

#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/spatial/force.hpp>
#include <pinocchio/spatial/motion.hpp>
#include <pinocchio/spatial/se3.hpp>

#include "galileo/common/meta/dimension.hpp"
#include "galileo/core/system-spec.hpp"

// Macros to import the types and constants from a multibody spec
#define GALILEO_MULTIBODY_SPEC_SCALARS_TYPEDEF(MultibodySpec) GALILEO_SYSTEM_SPEC_SCALARS_TYPEDEF(MultibodySpec::SS);

#define GALILEO_MULTIBODY_SPEC_CONSTANTS_TYPEDEF(MultibodySpec) \
    GALILEO_SYSTEM_SPEC_CONSTANTS_TYPEDEF(MultibodySpec::SS); \
    static constexpr int NQb = MultibodySpec::NQb; \
    static constexpr int NQj = MultibodySpec::NQj; \
    static constexpr int NVb = MultibodySpec::NVb; \
    static constexpr int NVj = MultibodySpec::NVj; \
    static constexpr int NRotors = MultibodySpec::NRotors; \
    using DimNQb_t = typename MultibodySpec::DimNQb_t; \
    using DimNQj_t = typename MultibodySpec::DimNQj_t; \
    using DimNVb_t = typename MultibodySpec::DimNVb_t; \
    using DimNVj_t = typename MultibodySpec::DimNVj_t; \
    using DimNRotors_t = typename MultibodySpec::DimNRotors_t;

#define GALILEO_MULTIBODY_SPEC_EIGEN_TYPES_TYPEDEF(MultibodySpec) \
    GALILEO_SYSTEM_SPEC_EIGEN_TYPES_TYPEDEF(MultibodySpec::SS); \
    using VectorNqb_t = typename MultibodySpec::VectorNqb_t; \
    using VectorNqj_t = typename MultibodySpec::VectorNqj_t; \
    using VectorNvb_t = typename MultibodySpec::VectorNvb_t; \
    using VectorNvj_t = typename MultibodySpec::VectorNvj_t;

#define GALILEO_MULTIBODY_SPECIFIC_TYPES_TYPEDEF(MultibodySpec) \
    using RobotModel_t = typename MultibodySpec::RobotModel_t; \
    using RobotData_t = typename MultibodySpec::RobotData_t; \
    using FrameIndex_t = typename MultibodySpec::FrameIndex_t; \
    using ReferenceFrame_t = typename MultibodySpec::ReferenceFrame_t; \
    using SE3_t = typename MultibodySpec::SE3_t; \
    using Motion_t = typename MultibodySpec::Motion_t; \
    using Force_t = typename MultibodySpec::Force_t; \
    using ActionMatrix_t = typename MultibodySpec::ActionMatrix_t;

#define GALILEO_MULTIBODY_SPEC_MASTER_TYPEDEF(MultibodySpec) \
    GALILEO_MULTIBODY_SPEC_SCALARS_TYPEDEF(MultibodySpec); \
    GALILEO_MULTIBODY_SPEC_CONSTANTS_TYPEDEF(MultibodySpec); \
    GALILEO_MULTIBODY_SPEC_EIGEN_TYPES_TYPEDEF(MultibodySpec); \
    GALILEO_MULTIBODY_SPECIFIC_TYPES_TYPEDEF(MultibodySpec);

namespace galileo
{

    /* ---------------------------------------------------------------- */
    /* Defines the multibody system dimensions and types. */
    /* ---------------------------------------------------------------- */
    template <typename BasicSpec,
              int _NQb,
              int _NQj,
              int _NVb,
              int _NVj,
              int _NRotors,
              template <typename> class StateTpl,
              template <typename> class ActuationTpl>
    struct MultibodySpecTpl
        : public SystemSpecTpl<BasicSpec, AddDim_v<_NQb, _NQj>, AddDim_v<_NVb, _NVj>, AddDim_v<_NVj, _NRotors>>
    {
        using MS = MultibodySpecTpl<BasicSpec, _NQb, _NQj, _NVb, _NVj, _NRotors, StateTpl, ActuationTpl>;
        using Base = SystemSpecTpl<BasicSpec, AddDim_v<_NQb, _NQj>, AddDim_v<_NVb, _NVj>, AddDim_v<_NVj, _NRotors>>;

        using BS = BasicSpec;
        GALILEO_SYSTEM_SPEC_MASTER_TYPEDEF(Base);

        /* ---------------------------------------------------------------- */
        /* Additional multibody dimensions */
        /* ---------------------------------------------------------------- */
        using DimNQb_t = DimensionTpl<_NQb>;         // Dimension of floating base generalized coordinates
        using DimNQj_t = DimensionTpl<_NQj>;         // Dimension of joint generalized coordinates
        using DimNVb_t = DimensionTpl<_NVb>;         // Dimension of floating base generalized velocities
        using DimNVj_t = DimensionTpl<_NVj>;         // Dimension of joint generalized velocities
        using DimNRotors_t = DimensionTpl<_NRotors>; // Number of rotors attached to the floating base

        static constexpr int NQb = DimNQb_t::Value;
        static constexpr int NQj = DimNQj_t::Value;
        static constexpr int NVb = DimNVb_t::Value;
        static constexpr int NVj = DimNVj_t::Value;
        static constexpr int NRotors = DimNRotors_t::Value;

        /* ---------------------------------------------------------------- */
        /* Essential Eigen types for multibody dimensions */
        /* ---------------------------------------------------------------- */
        using VectorNqb_t = Eigen::GMatrix<VarScalar, NQb, 1, Options>;
        using VectorNqj_t = Eigen::GMatrix<VarScalar, NQj, 1, Options>;
        using VectorNvb_t = Eigen::GMatrix<VarScalar, NVb, 1, Options>;
        using VectorNvj_t = Eigen::GMatrix<VarScalar, NVj, 1, Options>;

        /* ---------------------------------------------------------------- */
        /* Pinocchio-specific types */
        /* ---------------------------------------------------------------- */
        using RobotModel_t = pinocchio::ModelTpl<VarScalar, Options>;
        using RobotData_t = pinocchio::DataTpl<VarScalar, Options>;
        using FrameIndex_t = pinocchio::FrameIndex;
        using ReferenceFrame_t = pinocchio::ReferenceFrame;
        using SE3_t = pinocchio::SE3Tpl<VarScalar, Options>;
        using Motion_t = pinocchio::MotionTpl<VarScalar, Options>;
        using Force_t = pinocchio::ForceTpl<VarScalar, Options>;
        using ActionMatrix_t = typename SE3_t::ActionMatrixType;

        /* ---------------------------------------------------------------- */
        /* Template types */
        /* ---------------------------------------------------------------- */
        using State_t = StateTpl<MS>;

        using ActuationMeta_t = ActuationTpl<MS>;
        using ActuationModel_t = typename traits<ActuationMeta_t>::Model_t;
        using ActuationData_t = typename traits<ActuationMeta_t>::Data_t;

        /* ---------------------------------------------------------------- */
        /* Additional dimension storage */
        /* ---------------------------------------------------------------- */
        DimNQb_t nqb_dim_;
        DimNQj_t nqj_dim_;
        DimNVb_t nvb_dim_;
        DimNVj_t nvj_dim_;
        DimNRotors_t nrotors_dim_;

        using Base::nq_dim_;
        using Base::nv_dim_;
        using Base::nx_dim_;
        using Base::ndx_dim_;
        using Base::nua_dim_;

        /* ---------------------------------------------------------------- */
        /* Constructors */
        /* ---------------------------------------------------------------- */
        MultibodySpecTpl() : Base(), nqb_dim_{}, nqj_dim_{}, nvb_dim_{}, nvj_dim_{}, nrotors_dim_{} {}

        /* ---------------------------------------------------------------- */
        /* Accessors for multibody dimensions */
        /* ---------------------------------------------------------------- */
        const DimNQb_t &get_nqb_dim() const { return nqb_dim_; }
        int get_nqb() const { return nqb_dim_.value(); }
        const DimNQj_t &get_nqj_dim() const { return nqj_dim_; }
        int get_nqj() const { return nqj_dim_.value(); }
        const DimNVb_t &get_nvb_dim() const { return nvb_dim_; }
        int get_nvb() const { return nvb_dim_.value(); }
        const DimNVj_t &get_nvj_dim() const { return nvj_dim_; }
        int get_nvj() const { return nvj_dim_.value(); }
        const DimNRotors_t &get_nrotors_dim() const { return nrotors_dim_; }
        int get_nrotors() const { return nrotors_dim_.value(); }

        using Base::get_nq;
        using Base::get_nq_dim;
        using Base::get_nv;
        using Base::get_nv_dim;
        using Base::get_nx;
        using Base::get_nx_dim;
        using Base::get_ndx;
        using Base::get_ndx_dim;
        using Base::get_nua;
        using Base::get_nua_dim;

        inline bool is_valid_spec() const
        {
            bool valid_base = Base::is_valid_spec();
            bool valid_nq = (get_nq() == get_nqb() + get_nqj());
            bool valid_nv = (get_nv() == get_nvb() + get_nvj());
            bool valid_nua = (get_nua() == get_nvj() + get_nrotors());
            return valid_base && valid_nq && valid_nv && valid_nua;
        }

        void display(std::ostream &os, const std::string &indent = "  ") const
        {
            Base::display(os, indent);

            os << indent << "Multibody Dimensions: {\n";
            os << indent << "  Base Configuration:\n";
            os << indent << "    NQb (Base generalized coordinates):   " << get_nqb_dim() << "\n";
            os << indent << "    NVb (Base generalized velocities):    " << get_nvb_dim() << "\n";
            os << indent << "  Joint Configuration:\n";
            os << indent << "    NQj (Joint generalized coordinates):  " << get_nqj_dim() << "\n";
            os << indent << "    NVj (Joint generalized velocities):   " << get_nvj_dim() << "\n";
            os << indent << "  Additional Actuation:\n";
            os << indent << "    NRotors (Number of rotors):           " << get_nrotors_dim() << "\n";
            os << indent << "}\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const MultibodySpecTpl &ms)
        {
            os << "MultibodySpec: {\n";
            ms.display(os, "  ");
            os << "}";
            return os;
        }
    };

} // namespace galileo

#endif // __galileo_multibody_multibody_spec_hpp__
