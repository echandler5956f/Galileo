#ifndef __galileo_core_system_spec_hpp__
#define __galileo_core_system_spec_hpp__

#include "galileo/core/fwd.hpp"

// Macros to import the types and constants from a system spec
#define GALILEO_SYSTEM_SPEC_SCALARS_TYPEDEF(SystemSpec) GALILEO_BASIC_SPEC_SCALARS_TYPEDEF(SystemSpec::BS);

#define GALILEO_SYSTEM_SPEC_CONSTANTS_TYPEDEF(SystemSpec) \
    static constexpr int NQ = SystemSpec::NQ; \
    static constexpr int NV = SystemSpec::NV; \
    static constexpr int NX = SystemSpec::NX; \
    static constexpr int NDX = SystemSpec::NDX; \
    static constexpr int NUa = SystemSpec::NUa; \
    using DimNQ_t = typename SystemSpec::DimNQ_t; \
    using DimNV_t = typename SystemSpec::DimNV_t; \
    using DimNX_t = typename SystemSpec::DimNX_t; \
    using DimNDX_t = typename SystemSpec::DimNDX_t; \
    using DimNUa_t = typename SystemSpec::DimNUa_t;

#define GALILEO_SYSTEM_SPEC_EIGEN_TYPES_TYPEDEF(SystemSpec) \
    GALILEO_BASIC_SPEC_FIXED_SIZE_EIGEN_TYPES_TYPEDEF(SystemSpec::BS) \
    GALILEO_BASIC_SPEC_DYNAMIC_SIZE_EIGEN_TYPES_TYPEDEF(SystemSpec::BS) \
    using VectorNx_t = typename SystemSpec::VectorNx_t; \
    using VectorNua_t = typename SystemSpec::VectorNua_t; \
    using VectorNdx_t = typename SystemSpec::VectorNdx_t; \
    using VectorNq_t = typename SystemSpec::VectorNq_t; \
    using VectorNv_t = typename SystemSpec::VectorNv_t; \
    using MatrixNx_t = typename SystemSpec::MatrixNx_t; \
    using MatrixNua_t = typename SystemSpec::MatrixNua_t; \
    using MatrixNdx_t = typename SystemSpec::MatrixNdx_t; \
    using MatrixNq_t = typename SystemSpec::MatrixNq_t; \
    using MatrixNv_t = typename SystemSpec::MatrixNv_t; \
    using MatrixNvNdx_t = typename SystemSpec::MatrixNvNdx_t; \
    using MatrixNvNua_t = typename SystemSpec::MatrixNvNua_t; \
    using MatrixNuaNv_t = typename SystemSpec::MatrixNuaNv_t; \
    using MatrixNdxNua_t = typename SystemSpec::MatrixNdxNua_t; \
    using MatrixNuaNdx_t = typename SystemSpec::MatrixNuaNdx_t; \
    using MatrixNv6_t = typename SystemSpec::MatrixNv6_t; \
    using Matrix6Nv_t = typename SystemSpec::Matrix6Nv_t; \
    using Matrix6Ndx_t = typename SystemSpec::Matrix6Ndx_t;

#define GALILEO_SYSTEM_SPEC_MASTER_TYPEDEF(SystemSpec) \
    GALILEO_SYSTEM_SPEC_SCALARS_TYPEDEF(SystemSpec); \
    GALILEO_SYSTEM_SPEC_CONSTANTS_TYPEDEF(SystemSpec); \
    GALILEO_SYSTEM_SPEC_EIGEN_TYPES_TYPEDEF(SystemSpec);

namespace galileo
{

    /* ---------------------------------------------------------------- */
    /* Defines the system dimensions and types used in the core library. */
    /* ---------------------------------------------------------------- */
    template <typename BasicSpec, int _NQ, int _NV, int _NUa>
    struct SystemSpecTpl
    {
        using SS = SystemSpecTpl<BasicSpec, _NQ, _NV, _NUa>;

        /* ---------------------------------------------------------------- */
        /* Import the basic spec types and constants */
        /* ---------------------------------------------------------------- */
        using BS = BasicSpec;
        GALILEO_BASIC_SPEC_MASTER_TYPEDEF(BS);

        /* ---------------------------------------------------------------- */
        /* Dimensions and constants */
        /* ---------------------------------------------------------------- */
        using DimNQ_t = DimensionTpl<_NQ>;           // Dimension of generalized coordinates
        using DimNV_t = DimensionTpl<_NV>;           // Dimension of generalized velocities
        using DimNX_t = AddDim_t<DimNQ_t, DimNV_t>;  // Dimension of state
        using DimNDX_t = AddDim_t<DimNV_t, DimNV_t>; // Dimension of state tangent space
        using DimNUa_t = DimensionTpl<_NUa>;         // Dimension of actuated inputs

        static constexpr int NQ = DimNQ_t::Value;
        static constexpr int NV = DimNV_t::Value;
        static constexpr int NX = DimNX_t::Value;
        static constexpr int NDX = DimNDX_t::Value;
        static constexpr int NUa = DimNUa_t::Value;

        /* ---------------------------------------------------------------- */
        /* Essential Eigen types for system dimensions */
        /* ---------------------------------------------------------------- */
        using VectorNx_t = Eigen::GMatrix<VarScalar, NX, 1, Options>;
        using VectorNua_t = Eigen::GMatrix<VarScalar, NUa, 1, Options>;
        using VectorNdx_t = Eigen::GMatrix<VarScalar, NDX, 1, Options>;
        using VectorNq_t = Eigen::GMatrix<VarScalar, NQ, 1, Options>;
        using VectorNv_t = Eigen::GMatrix<VarScalar, NV, 1, Options>;

        using MatrixNx_t = Eigen::GMatrix<VarScalar, NX, NX, Options>;
        using MatrixNua_t = Eigen::GMatrix<VarScalar, NUa, NUa, Options>;
        using MatrixNdx_t = Eigen::GMatrix<VarScalar, NDX, NDX, Options>;
        using MatrixNq_t = Eigen::GMatrix<VarScalar, NQ, NQ, Options>;
        using MatrixNv_t = Eigen::GMatrix<VarScalar, NV, NV, Options>;

        using MatrixNvNdx_t = Eigen::GMatrix<VarScalar, NV, NDX, Options>;
        using MatrixNvNua_t = Eigen::GMatrix<VarScalar, NV, NUa, Options>;
        using MatrixNuaNv_t = Eigen::GMatrix<VarScalar, NUa, NV, Options>;
        using MatrixNdxNua_t = Eigen::GMatrix<VarScalar, NDX, NUa, Options>;
        using MatrixNuaNdx_t = Eigen::GMatrix<VarScalar, NUa, NDX, Options>;

        using MatrixNv6_t = Eigen::GMatrix<VarScalar, NV, 6, Options>;
        using Matrix6Nv_t = Eigen::GMatrix<VarScalar, 6, NV, Options>;
        using Matrix6Ndx_t = Eigen::GMatrix<VarScalar, 6, NDX, Options>;

        /* ---------------------------------------------------------------- */
        /* Dimension storage */
        /* ---------------------------------------------------------------- */
        DimNQ_t nq_dim_;
        DimNV_t nv_dim_;
        DimNX_t nx_dim_;
        DimNDX_t ndx_dim_;
        DimNUa_t nua_dim_;

        /* ---------------------------------------------------------------- */
        /* Accessors for the SystemSpec dimensions */
        /* ---------------------------------------------------------------- */
        const DimNQ_t &get_nq_dim() const { return nq_dim_; }
        int get_nq() const { return nq_dim_.value(); }
        const DimNV_t &get_nv_dim() const { return nv_dim_; }
        int get_nv() const { return nv_dim_.value(); }
        const DimNX_t &get_nx_dim() const { return nx_dim_; }
        int get_nx() const { return nx_dim_.value(); }
        const DimNDX_t &get_ndx_dim() const { return ndx_dim_; }
        int get_ndx() const { return ndx_dim_.value(); }
        const DimNUa_t &get_nua_dim() const { return nua_dim_; }
        int get_nua() const { return nua_dim_.value(); }

        inline bool is_valid_spec() const
        {
            // Basic consistency checks
            bool valid_nx = (get_nx() == get_nq() + get_nv());
            bool valid_ndx = (get_ndx() == 2 * get_nv());
            return valid_nx && valid_ndx;
        }

        void display(std::ostream &os, const std::string &indent = "  ") const
        {
            os << indent << "BasicSpec: {";
            BS{}.display(os, "");
            os << "}\n";

            os << indent << "System Dimensions: {\n";
            os << indent << "  Configuration:\n";
            os << indent << "    NQ (Generalized coordinates):  " << get_nq_dim() << "\n";
            os << indent << "    NV (Generalized velocities):   " << get_nv_dim() << "\n";
            os << indent << "  Actuation:\n";
            os << indent << "    NUa (Actuated inputs):         " << get_nua_dim() << "\n";
            os << indent << "  State:\n";
            os << indent << "    NX (State dimension):          " << get_nx_dim() << "\n";
            os << indent << "    NDX (State tangent dimension): " << get_ndx_dim() << "\n";
            os << indent << "}\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const SystemSpecTpl &ss)
        {
            os << "SystemSpec: {\n";
            ss.display(os, "  ");
            os << "}";
            return os;
        }

    protected:
        inline SystemSpecTpl(int NQ_, int NV_, int NUa_)
            : nq_dim_{NQ_}, nv_dim_{NV_}, nx_dim_{NQ_ + NV_}, ndx_dim_{2 * NV_}, nua_dim_{NUa_}
        {
        }
        inline SystemSpecTpl(const SystemSpecTpl &clone)
            : nq_dim_{clone.nq_dim_},
              nv_dim_{clone.nv_dim_},
              nx_dim_{clone.nx_dim_},
              ndx_dim_{clone.ndx_dim_},
              nua_dim_{clone.nua_dim_}
        {
        }
        inline SystemSpecTpl &operator=(const SystemSpecTpl &clone)
        {
            nq_dim_ = clone.nq_dim_;
            nv_dim_ = clone.nv_dim_;
            nx_dim_ = clone.nx_dim_;
            ndx_dim_ = clone.ndx_dim_;
            nua_dim_ = clone.nua_dim_;
            return *this;
        }
    };

    /* ---------------------------------------------------------------- */
    /* Defines the default spec based on State and Actuation types. */
    /* ---------------------------------------------------------------- */
    template <typename BasicSpec,
              int _NQ,
              int _NV,
              int _NUa,
              template <typename> class StateTpl,
              template <typename> class ActuationTpl>
    struct DefaultSpecTpl : public SystemSpecTpl<BasicSpec, _NQ, _NV, _NUa>
    {
        using DS = DefaultSpecTpl<BasicSpec, _NQ, _NV, _NUa, StateTpl, ActuationTpl>;
        using SS = SystemSpecTpl<BasicSpec, _NQ, _NV, _NUa>;
        using Base = SS;
        using BS = BasicSpec;

        /* ---------------------------------------------------------------- */
        /* Import the system spec types and constants */
        /* ---------------------------------------------------------------- */
        GALILEO_SYSTEM_SPEC_MASTER_TYPEDEF(Base);

        /* ---------------------------------------------------------------- */
        /* Template types */
        /* ---------------------------------------------------------------- */
        using State_t = StateTpl<DS>;

        using ActuationMeta_t = ActuationTpl<DS>;
        using ActuationModel_t = typename traits<ActuationMeta_t>::Model_t;
        using ActuationData_t = typename traits<ActuationMeta_t>::Data_t;

        using Base::nq_dim_;
        using Base::nv_dim_;
        using Base::nx_dim_;
        using Base::ndx_dim_;
        using Base::nua_dim_;

        DefaultSpecTpl(int NQ_, int NV_, int NUa_) : Base(NQ_, NV_, NUa_) {}

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
        using Base::is_valid_spec;
        using Base::display;

        friend std::ostream &operator<<(std::ostream &os, const DefaultSpecTpl &ds)
        {
            os << "DefaultSpec: {\n";
            ds.display(os, "  ");
            os << "}";
            return os;
        }
    };

} // namespace galileo

#endif // __galileo_core_system_spec_hpp__
