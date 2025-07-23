#ifndef __galileo_predictive_nodes_node_contact_fwddyn_hpp__
#define __galileo_predictive_nodes_node_contact_fwddyn_hpp__

#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/contact-dynamics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <pinocchio/algorithm/rnea.hpp>

#include "galileo/common/math/matrix-decomposition.hpp"

#include "galileo/predictive/nodes/node-base.hpp"

#include "galileo/multibody/contacts/contact-manager.hpp"
#include "galileo/multibody/contacts/fwd.hpp"

#include "galileo/core/data/data-collector-default.hpp"

#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct NodeContactFwdDynTpl;

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<NodeContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = NodeModelContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Data_t = NodeDataContactFwdDynTpl<PS, ContactCollectionTpl>;

        using DimNC_t = DimensionTpl<Eigen::Dynamic>;
        using DimNU_t = typename PS::DimNUa_t;

        using ContactManagerMeta_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using ContactModelManager_t = typename traits<ContactManagerMeta_t>::ModelManager_t;
        using ContactDataManager_t = typename traits<ContactManagerMeta_t>::DataManager_t;

        using DataCollector_t = DataCollectorDefaultTpl<PS, ContactCollectionTpl>;

        // using NC = traits<NodeDerived>::NC;
        // using Kinv_t = Eigen::GMatrix<VarScalar, PS::NV + NC, PS::NV + NC>;
        // using Jstatic_t = Eigen::GMatrix<VarScalar, PS::NV, NU + NC>;
        using Kinv_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, Eigen::Dynamic>;
        using MatrixNcNdx_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNDX_t::Value>;
        using MatrixNcNu_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, DimNU_t::Value>;
        using Jstatic_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNV_t::Value, Eigen::Dynamic>;

        using MatrixNvNc_t = Eigen::GMatrix<typename PS::VarScalar, PS::DimNV_t::Value, Eigen::Dynamic>;
        using MatrixNcNv_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNV_t::Value>;
        using MatrixNc_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, Eigen::Dynamic>;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<NodeDataContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct traits<NodeModelContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct NodeDataContactFwdDynTpl
        : public NodeDataBase<NodeDataContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = NodeDataBase<NodeDataContactFwdDynTpl<PS, ContactCollectionTpl>, PS>;

        using ContactManagerMeta_t = typename traits<Meta_t>::ContactManagerMeta_t;
        using ContactModelManager_t = typename traits<Meta_t>::ContactModelManager_t;
        using ContactDataManager_t = typename traits<Meta_t>::ContactDataManager_t;

        using DataCollector_t = typename traits<Meta_t>::DataCollector_t;
        using JointData_t = typename DataCollector_t::JointData_t;

        using Kinv_t = typename traits<Meta_t>::Kinv_t;
        using MatrixNcNdx_t = typename traits<Meta_t>::MatrixNcNdx_t;
        using MatrixNcNu_t = typename traits<Meta_t>::MatrixNcNu_t;
        using Jstatic_t = typename traits<Meta_t>::Jstatic_t;

        DEFAULT_ACCESSOR(CostDataManager_t, costs);
        DEFAULT_ACCESSOR(ConstraintDataManager_t, constraints);
        DEFAULT_ACCESSOR(ActuationData_t, actuation);

        DEFAULT_ACCESSOR(XAcc_t, XAcc);
        DEFAULT_ACCESSOR(XAccx_t, XAccx);
        DEFAULT_ACCESSOR(XAccu_t, XAccu);

        // Notice that we are overriding the default accessors for the costs and constraints

        L_t &L_accessor() { return costs.L; }
        const L_t &L_accessor() const { return costs.L; }

        Lx_t &Lx_accessor() { return costs.Lx; }
        const Lx_t &Lx_accessor() const { return costs.Lx; }

        Lu_t &Lu_accessor() { return costs.Lu; }
        const Lu_t &Lu_accessor() const { return costs.Lu; }

        Lxx_t &Lxx_accessor() { return costs.Lxx; }
        const Lxx_t &Lxx_accessor() const { return costs.Lxx; }

        Lxu_t &Lxu_accessor() { return costs.Lxu; }
        const Lxu_t &Lxu_accessor() const { return costs.Lxu; }

        Luu_t &Luu_accessor() { return costs.Luu; }
        const Luu_t &Luu_accessor() const { return costs.Luu; }

        H_t &H_accessor() { return constraints.H; }
        const H_t &H_accessor() const { return constraints.H; }

        Hx_t &Hx_accessor() { return constraints.Hx; }
        const Hx_t &Hx_accessor() const { return constraints.Hx; }

        Hu_t &Hu_accessor() { return constraints.Hu; }
        const Hu_t &Hu_accessor() const { return constraints.Hu; }

        // For now the inequality constraints accessor points to the equality constraints,
        // until we implement inequality constraints
        G_t &G_accessor() { return constraints.H; }
        const G_t &G_accessor() const { return constraints.H; }

        Gx_t &Gx_accessor() { return constraints.Hx; }
        const Gx_t &Gx_accessor() const { return constraints.Hx; }

        Gu_t &Gu_accessor() { return constraints.Hu; }
        const Gu_t &Gu_accessor() const { return constraints.Hu; }

        NodeDataContactFwdDynTpl(const Model_t &model)
            : XAcc(model.get_ps().get_nv()),
              XAccx(model.get_ps().get_nv(), model.get_ps().get_ndx()),
              XAccu(model.get_ps().get_nv(), model.get_ps().get_nu()),
              robot(RobotData_t(model.get_robot())),
              actuation(model.get_actuation().createData()),
              joint(JointData_t(model.get_ps())),
              contacts(model.get_contacts().createData(&robot)),
              data_collector(&robot, &actuation, &joint, &contacts),
              costs(model.get_costs().createData(&data_collector)),
              constraints(model.get_constraints().createData(&data_collector)),
              Kinv(model.get_ps().get_nv() +
                       model.get_contacts().get_nc_total(),
                   model.get_ps().get_nv() +
                       model.get_contacts().get_nc_total()),
              df_dx(model.get_contacts().get_nc_total(),
                    model.get_ps().get_ndx()),
              df_du(model.get_contacts().get_nc_total(), model.get_ps().get_nu()),
              tmp_xstatic(model.get_ps().get_nx()),
              tmp_Jstatic(model.get_ps().get_nv(),
                          model.get_ps().get_nu() + model.get_contacts().get_nc_total())
        {
            XAcc.setZero();
            XAccx.setZero();
            XAccu.setZero();
            joint.dtau_du.diagonal().setOnes();
            Kinv.setZero();
            df_dx.setZero();
            df_du.setZero();
            tmp_xstatic.setZero();
            tmp_Jstatic.setZero();
            robot.lambda_c.resize(model.get_contacts().get_nc_total());
            robot.lambda_c.setZero();
        }

        XAcc_t XAcc;
        XAccx_t XAccx;
        XAccu_t XAccu;

        RobotData_t robot;
        ActuationData_t actuation;
        JointData_t joint;
        ContactDataManager_t contacts;
        DataCollector_t data_collector;
        CostDataManager_t costs;
        ConstraintDataManager_t constraints;

        Kinv_t Kinv;
        MatrixNcNdx_t df_dx;
        MatrixNcNu_t df_du;
        VectorNx_t tmp_xstatic;
        Jstatic_t tmp_Jstatic;

    }; // class NodeDataContactFwdDynTpl

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    class NodeModelContactFwdDynTpl
        : public NodeModelBase<NodeModelContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = NodeModelBase<NodeModelContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>;

        using ContactManagerMeta_t = typename traits<Meta_t>::ContactManagerMeta_t;
        using ContactModelManager_t = typename traits<Meta_t>::ContactModelManager_t;
        using ContactDataManager_t = typename traits<Meta_t>::ContactDataManager_t;

        using MatrixNvNc_t = typename traits<Meta_t>::MatrixNvNc_t;
        using MatrixNcNv_t = typename traits<Meta_t>::MatrixNcNv_t;
        using MatrixNc_t = typename traits<Meta_t>::MatrixNc_t;

        NodeModelContactFwdDynTpl(PS &ps, const CostModelManager_t &costs,
                                  const ConstraintModelManager_t &constraints,
                                  const ContactModelManager_t &contacts,
                                  const ActuationModel_t &actuation,
                                  const NumScalar &JMinvJt_damping,
                                  const bool enable_force)
            : Base(ps),
              costs_(costs),
              constraints_(constraints),
              contacts_(contacts),
              actuation_(actuation),
              with_armature_(true),
              armature_(get_ps().get_nv()),
              JMinvJt_damping_(fabs(JMinvJt_damping)),
              enable_force_(enable_force)
        {
            int nua = get_ps().get_nua();
            get_ps().set_nu(nua);
            set_u_lb(NumScalar(-1.) * get_robot().effortLimit.tail(nua));
            set_u_ub(NumScalar(+1.) * get_robot().effortLimit.tail(nua));

            armature_.setZero();
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            auto nc_active_dim = get_contacts().get_nc_active_dim();

            const auto q = head(x, get_ps().get_nq_dim());
            const auto v = tail(x, get_ps().get_nv_dim());
            pinocchio::computeAllTerms(get_robot(), data.robot, q, v);
            pinocchio::computeCentroidalMomentum(get_robot(), data.robot);

            if (!with_armature_)
            {
                data.robot.M.diagonal() += armature_;
            }
            get_actuation().calc(data.actuation, x, u);
            get_contacts().calc(data.contacts, x);
            pinocchio::forwardDynamics(
                get_robot(), data.robot, data.actuation.tau,
                topRows(data.contacts.Jc, nc_active_dim),
                head(data.contacts.a0, nc_active_dim),
                JMinvJt_damping_);
            data.XAcc = data.robot.ddq;
            get_contacts().updateAcceleration(data.contacts, data.robot.ddq);
            get_contacts().updateForce(data.contacts, data.robot.lambda_c);
            data.joint.a = data.robot.ddq;
            data.joint.tau = u;
            get_costs().calc(data.costs, x, u);
            data.L_accessor() = data.costs.L;
            if (get_constraints().get_nh() > 0 || get_constraints().get_nh() > 0)
            {
                // data.constraints.resize(this, data);
                get_constraints().calc(data.constraints, x, u);
            }
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const auto q = head(x, get_ps().get_nq_dim());
            const auto v = tail(x, get_ps().get_nv_dim());
            pinocchio::computeAllTerms(get_robot(), data.robot, q, v);
            pinocchio::computeCentroidalMomentum(get_robot(), data.robot);
            get_costs().calc(data.costs, x);
            data.L_accessor() = data.costs.L;
            if (get_constraints().get_nh() > 0 || get_constraints().get_nh() > 0)
            {
                // data.constraints.resize(this, data);
                get_constraints().calc(data.constraints, x);
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            auto nc_active_dim = get_contacts().get_nc_active_dim();
            const auto q = head(x, get_ps().get_nq_dim());
            const auto v = tail(x, get_ps().get_nv_dim());

            // Computing the dynamics derivatives
            // We resize the Kinv matrix because Eigen cannot call block operations
            // recursively: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=408. Therefore,
            // it is not possible to pass data.Kinv.topLeftCorner(nv + nc, nv + nc)
            data.Kinv.resize(get_ps().get_nv() + nc_active_dim.value(), get_ps().get_nv() + nc_active_dim.value());
            pinocchio::computeRNEADerivatives(get_robot(), data.robot, q, v, data.XAcc,
                                              data.contacts.fext);
            get_contacts().updateRneaDiff(data.contacts, data.robot);
            pinocchio::getKKTContactDynamicMatrixInverse(
                get_robot(), data.robot, topRows(data.contacts.Jc, nc_active_dim), data.Kinv);

            get_actuation().calcDiff(data.actuation, x, u);
            get_contacts().calcDiff(data.contacts, x);

            const auto a_partial_dtau = topLeftCorner(data.Kinv, get_ps().get_nv_dim(), get_ps().get_nv_dim());
            const auto a_partial_da = topRightCorner(data.Kinv, get_ps().get_nv_dim(), nc_active_dim);
            const auto f_partial_dtau = bottomLeftCorner(data.Kinv, nc_active_dim, get_ps().get_nv_dim());
            const auto f_partial_da = bottomRightCorner(data.Kinv, nc_active_dim, nc_active_dim);

            leftCols(data.XAccx, get_ps().get_nv_dim()).noalias() = -a_partial_dtau * data.robot.dtau_dq;
            rightCols(data.XAccx, get_ps().get_nv_dim()).noalias() = -a_partial_dtau * data.robot.dtau_dv;
            data.XAccx.noalias() -= a_partial_da * topRows(data.contacts.da0_dx, nc_active_dim);
            data.XAccx.noalias() += a_partial_dtau * data.actuation.dtau_dx;
            data.XAccu.noalias() = a_partial_dtau * data.actuation.dtau_du;
            data.joint.da_dx = data.XAccx;
            data.joint.da_du = data.XAccu;

            // Computing the cost derivatives
            if (enable_force_)
            {
                topLeftCorner(data.df_dx, nc_active_dim, get_ps().get_nv_dim()).noalias() =
                    f_partial_dtau * data.robot.dtau_dq;
                topRightCorner(data.df_dx, nc_active_dim, get_ps().get_nv_dim()).noalias() =
                    f_partial_dtau * data.robot.dtau_dv;
                topRows(data.df_dx, nc_active_dim).noalias() +=
                    f_partial_da * topRows(data.contacts.da0_dx, nc_active_dim);
                topRows(data.df_dx, nc_active_dim).noalias() -=
                    f_partial_dtau * data.actuation.dtau_dx;
                topRows(data.df_du, nc_active_dim).noalias() =
                    -f_partial_dtau * data.actuation.dtau_du;
                get_contacts().updateAccelerationDiff(data.contacts,
                                                      bottomRows(data.XAccx, get_ps().get_nv_dim()));
                get_contacts().updateForceDiff(data.contacts, topRows(data.df_dx, nc_active_dim),
                                               topRows(data.df_du, nc_active_dim));
            }
            get_costs().calcDiff(data.costs, x, u);
            if (get_constraints().get_nh() > 0 || get_constraints().get_nh() > 0)
            {
                get_constraints().calcDiff(data.constraints, x, u);
            }
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            get_costs().calcDiff(data.costs, x);
            if (get_constraints().get_nh() > 0 || get_constraints().get_nh() > 0)
            {
                get_constraints().calcDiff(data.constraints, x);
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void quasiStatic(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlVectorType> &u,
                         const int maxiter, const NumScalar tol) const
        {
            const auto q = head(x, get_ps().get_nq_dim());
            auto nc_active_dim = get_contacts().get_nc_active_dim();

            head(data.tmp_xstatic, get_ps().get_nq_dim()) = q;
            tail(data.tmp_xstatic, get_ps().get_nv_dim()).setZero();
            u.setZero();

            pinocchio::computeAllTerms(get_robot(), data.robot, q,
                                       tail(data.tmp_xstatic, get_ps().get_nv_dim()));
            pinocchio::computeJointJacobians(get_robot(), data.robot, q);
            pinocchio::rnea(get_robot(), data.robot, q,
                            tail(data.tmp_xstatic, get_ps().get_nv_dim()),
                            tail(data.tmp_xstatic, get_ps().get_nv_dim()));
            get_actuation().calc(data.actuation, data.tmp_xstatic, u);
            get_actuation().calcDiff(data.actuation, data.tmp_xstatic, u);
            get_contacts().calc(data.contacts, data.tmp_xstatic);

            // Allocates memory
            data.tmp_Jstatic.conservativeResize(get_ps().get_nv(), get_ps().get_nu() + nc_active_dim.value());
            data.tmp_Jstatic.leftCols(get_ps().get_nu()) = data.actuation.dtau_du;
            rightCols(data.tmp_Jstatic, nc_active_dim) =
                topRows(data.contacts.Jc, nc_active_dim).transpose();

            u.noalias() = head(pseudoInverse(data.tmp_Jstatic) * data.robot.tau, get_ps().get_nu_dim());
            data.robot.tau.setZero();
        }

        Data_t createData() const
        {
            return Data_t(*this);
        }

        const CostModelManager_t &get_costs() const
        {
            return costs_.get();
        }

        const ConstraintModelManager_t &get_constraints() const
        {
            return constraints_.get();
        }

        const ContactModelManager_t &get_contacts() const
        {
            return contacts_.get();
        }

        const ActuationModel_t &get_actuation() const
        {
            return actuation_.get();
        }

        using Base::get_ps;

        using Base::get_robot;
        using Base::get_state;

        using Base::get_u_lb;
        using Base::get_u_ub;

        using Base::set_u_lb;
        using Base::set_u_ub;

    protected:
        std::reference_wrapper<const CostModelManager_t> costs_;
        std::reference_wrapper<const ConstraintModelManager_t> constraints_;
        std::reference_wrapper<const ContactModelManager_t> contacts_;
        std::reference_wrapper<const ActuationModel_t> actuation_;

        bool with_armature_ = true;
        VectorNv_t armature_;
        NumScalar JMinvJt_damping_;
        bool enable_force_ = false;

    }; // class NodeModelContactFwdDynTpl

} // namespace galileo

#endif // __galileo_predictive_nodes_node_contact_fwddyn_hpp__
