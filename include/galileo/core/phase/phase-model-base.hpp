#ifndef __galileo_core_phase_model_base_hpp__
#define __galileo_core_phase_model_base_hpp__

#include "galileo/core/phase/phase-base.hpp"
#include "galileo/utils/aligned-vector.hpp"
#include <limits>

#define GALILEO_PHASE_MODEL_TYPEDEF_GENERIC(Phase, TYPENAME)               \
    typedef TYPENAME traits<Phase>::Scalar Scalar;                           \
    typedef TYPENAME traits<Phase>::PhaseModelDerived PhaseModelDerived; \
    typedef TYPENAME traits<Phase>::PhaseDataDerived PhaseDataDerived;   \
    typedef TYPENAME traits<Phase>::NodeModel_t NodeModel_t;                 \
    typedef TYPENAME traits<Phase>::NodeData_t NodeData_t;                   \
    typedef TYPENAME traits<Phase>::State_t State_t;                         \
    typedef TYPENAME traits<Phase>::Control_t Control_t;                     \
    enum                                                                       \
    {                                                                          \
        Options = traits<Phase>::Options,                                    \
        NX = traits<Phase>::NX,                                              \
        NDX = traits<Phase>::NDX,                                            \
        NU = traits<Phase>::NU,                                              \
        NH = traits<Phase>::NH,                                              \
        NG = traits<Phase>::NG,                                              \
        NC = traits<Phase>::NC                                               \
    };                                                                         \
    typedef TYPENAME traits<Phase>::X_t X_t;                                 \
    typedef TYPENAME traits<Phase>::dX_t dX_t;                               \
    typedef TYPENAME traits<Phase>::U_t U_t;                                 \
    typedef TYPENAME traits<Phase>::H_t H_t;                                 \
    typedef TYPENAME traits<Phase>::G_t G_t;                                 \
    typedef TYPENAME traits<Phase>::C_t C_t;

#define GALILEO_PHASE_TYPEDEF_TEMPLATE(Phase) \
    GALILEO_PHASE_MODEL_TYPEDEF_GENERIC(Phase, typename)

#define GALILEO_PHASE_CAST_TYPE_SPECIALIZATION(PhaseModelTpl) \
    template <typename Scalar, typename NewScalar>                \
    struct CastType<NewScalar, PhaseModelTpl<Scalar>>           \
    {                                                             \
        typedef PhaseModelTpl<NewScalar> type;                  \
    }

namespace galileo
{
    
    template <typename Derived>
    class PhaseModelBase : CRTP<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PhaseDerived = typename traits<Derived>::PhaseDerived;
        GALILEO_PHASE_TYPEDEF_TEMPLATE(PhaseDerived);

        PhaseDataDerived createData() const
        {
            return derived().createData();
        }

        template <typename StateMatrixType, typename ControlMatrixType>
        void calc(PhaseDataDerived &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs,
                  const Eigen::MatrixBase<ControlMatrixType> &us) const
        {
            derived().calc(data, xs.derived(), us.derived());
        }

        template <typename StateMatrixType, typename ControlMatrixType>
        void calc(PhaseDataDerived &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs) const
        {
            derived().calc(data, xs.derived());
        }

        template <typename StateMatrixType, typename ControlMatrixType>
        void calcDiff(PhaseDataDerived &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlMatrixType> &us) const
        {
            derived().calcDiff(data, xs.derived(), us.derived());
        }

        template <typename StateMatrixType, typename ControlMatrixType>
        void calcDiff(PhaseDataDerived &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs) const
        {
            derived().calcDiff(data, xs.derived());
        }

        template <typename NewScalar>
        typename CastType<NewScalar, Derived>::type cast() const
        {
            return derived().template cast<NewScalar>();
        }

    protected:
        // Default constructor: protected.
        // Prevent the construction of stand-alone PhaseModelBase.
        inline PhaseModelBase() : period_(std::numeric_limits<double>::quiet_NaN), num_nodes_(std::numeric_limits<Eigen::Index>::max()), integ_constraint_size_(std::numeric_limits<Eigen::Index>::max())
        {
        }

        // Copy constructor: protected.
        // Copy of stand-alone PhaseModelBase are prevented, but can be used from inheriting
        // objects. Copy is done by calling copy operator.
        inline PhaseModelBase(const PhaseModelBase &clone)
        {
            *this = clone;
        }

        // Copy operator: protected.
        // Copy of stand-alone PhaseModelBase are prevented, but can be used from inheriting
        // objects.
        inline PhaseModelBase &operator=(const PhaseModelBase &clone)
        {
            node_models_ = clone.node_models_;
            state_desc_ = clone.state_desc_;
            control_desc_ = clone.control_desc_;
            period_ = clone.period_;
            num_nodes_ = clone.num_nodes_;
            nc_ = clone.nc_;
            return *this;
        }

        GALILEO_ALIGNED_STD_VECTOR(NodeModel_t)
        node_models_; // node models

        State_t state_desc_;     // state description
        Control_t control_desc_; // control description

        Scalar period_;    // time period
        Eigen::Index num_nodes_; // number of nodes
        Eigen::Index nc_;        // number of integration constraints

    }; // class PhaseModelBase

} // namespace galileo

#endif // __galileo_core_phase_model_base_hpp__
