#ifndef __galileo_predictive_nodes_node_contact_fwddyn_hpp__
#define __galileo_predictive_nodes_node_contact_fwddyn_hpp__

#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/contact-dynamics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <pinocchio/algorithm/rnea.hpp>

#include "galileo/predictive/nodes/node-base.hpp"

#include "galileo/multibody/contacts/contact-manager.hpp"
#include "galileo/multibody/contacts/fwd.hpp"

#include "galileo/core/data/data-collector-default.hpp"

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
        using MatrixNcNu_t = Eigen::GMatrix<typename PS::VarScalar, Eigen::Dynamic, PS::DimNU_t::Value>;
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
        using NodeData_t = typename traits<Meta_t>::NodeData_t;
        using NodeModel_t = typename traits<Meta_t>::NodeModel_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    struct NodeDataContactFwdDynTpl
        : public NodeDataBase<NodeDataContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using ContactManagerMeta_t = typename traits<Meta_t>::ContactManagerMeta_t;
        using ContactModelManager_t = typename traits<Meta_t>::ContactModelManager_t;
        using ContactDataManager_t = typename traits<Meta_t>::ContactDataManager_t;

        using DataCollector_t = typename traits<Meta_t>::DataCollector_t;
        using JointData_t = typename DataCollector_t::JointData_t;

        using Kinv_t = typename traits<Meta_t>::Kinv_t;
        using MatrixNcNdx_t = typename traits<Meta_t>::MatrixNcNdx_t;
        using MatrixNcNu_t = typename traits<Meta_t>::MatrixNcNu_t;
        using Jstatic_t = typename traits<Meta_t>::Jstatic_t;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

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

        template <typename DataCollector>
        NodeDataContactFwdDynTpl(const Model_t &model)
            : XAcc(model.get_ps().get_nv()),
              XAccx(model.get_ps().get_nv(), model.get_ps().get_ndx()),
              XAccu(model.get_ps().get_nv(), model.get_ps().get_nu()),
              robot(RobotData_t(model.get_robot())),
              actuation(model.get_actuation().createData()),
              joint(JointData_t(model.get_ps())),
              contacts(model.get_contacts().createData(&robot)),
              data_collector(&robot, &actuation, &joint, &contacts),
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
            constraints = model.get_constraints().createData(&data_collector);
            costs = model.get_costs().createData(&data_collector);
            Kinv.setZero();
            df_dx.setZero();
            df_du.setZero();
            tmp_xstatic.setZero();
            tmp_Jstatic.setZero();
            robot.lambda_c.resize(model.get_contacts().get_nc_total());
            robot.lambda_c.setZero();
        }

        CostDataManager_t costs;
        ConstraintDataManager_t constraints;
        ContactDataManager_t contacts;
        ActuationData_t actuation;
        JointData_t joint;
        RobotData_t robot;

        DataCollector_t data_collector;

        XAcc_t XAcc;
        XAccx_t XAccx;
        XAccu_t XAccu;

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
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            auto nc_active_dim = contacts_.get_nc_active_dim();

            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, PS::DimNQ_t::Value> q =
                head(x, get_ps().get_nq_dim());
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, PS::DimNV_t::Value> v =
                tail(x, get_ps().get_nv_dim());

            // Computing the forward dynamics with the holonomic constraints defined by
            // the contact model
            pinocchio::computeAllTerms(get_robot(), data.robot, q, v);
            pinocchio::computeCentroidalMomentum(get_robot(), data.robot);

            if (!with_armature_)
            {
                data.robot.M.diagonal() += armature_;
            }
            actuation_.calc(data.actuation, x, u);
            contacts_.calc(data.contacts, x);

            pinocchio::forwardDynamics(
                get_robot(), data.robot, data.actuation.tau,
                topRows(data.contacts.Jc, nc_active_dim),
                head(data.contacts.a0, nc_active_dim),
                JMinvJt_damping_);
            data.XAcc = data.robot.ddq;
            contacts_.updateAcceleration(data.contacts, data.robot.ddq);
            contacts_.updateForce(data.contacts, data.robot.lambda_c);
            data.joint.a = data.robot.ddq;
            data.joint.tau = u;
            costs_.calc(data.costs, x, u);
            data.cost = data.costs.cost;
            if (constraints_.get_ng() > 0 || constraints_.get_nh() > 0)
            {
                data.constraints.resize(this, data);
                constraints_.calc(data.constraints, x, u);
            }
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, PS::DimNQ_t::Value> q =
                head(x, get_ps().get_nq_dim());
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, PS::DimNV_t::Value> v =
                tail(x, get_ps().get_nv_dim());

            pinocchio::computeAllTerms(get_robot(), data.robot, q, v);
            pinocchio::computeCentroidalMomentum(get_robot(), data.robot);
            costs_.calc(data.costs, x);
            data.cost = data.costs.cost;
            if (constraints_.get_ng() > 0 || constraints_.get_nh() > 0)
            {
                data.constraints.resize(this, data);
                constraints_.calc(data.constraints, x);
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            auto nc_active_dim = contacts_.get_nc_active_dim();
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, PS::DimNQ_t::Value> q =
                head(x, get_ps().get_nq_dim());
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, PS::DimNV_t::Value> v =
                tail(x, get_ps().get_nv_dim());

            // Computing the dynamics derivatives
            // We resize the Kinv matrix because Eigen cannot call block operations
            // recursively: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=408. Therefore,
            // it is not possible to pass data.Kinv.topLeftCorner(nv + nc, nv + nc)
            data.Kinv.resize(get_ps().get_nv() + nc_active_dim.value(), get_ps().get_nv() + nc_active_dim.value());
            pinocchio::computeRNEADerivatives(get_robot(), data.robot, q, v, data.XAcc,
                                              data.contacts.fext);
            contacts_.updateRneaDiff(data.contacts, data.robot);
            pinocchio::getKKTContactDynamicMatrixInverse(
                get_robot(), data.robot, topRows(data.contacts.Jc, nc_active_dim), data.Kinv);

            actuation_.calcDiff(data.actuation, x, u);
            contacts_.calcDiff(data.contacts, x);

            const Eigen::Block<MatrixNv_t> a_partial_dtau = topLeftCorner(data.Kinv, get_ps().get_nv_dim(), get_ps().get_nv_dim());
            const Eigen::Block<MatrixNvNc_t> a_partial_da = topRightCorner(data.Kinv, get_ps().get_nv_dim(), nc_active_dim);
            const Eigen::Block<MatrixNcNv_t> f_partial_dtau = bottomLeftCorner(data.Kinv, nc_active_dim, get_ps().get_nv_dim());
            const Eigen::Block<MatrixNc_t> f_partial_da = bottomRightCorner(data.Kinv, nc_active_dim, nc_active_dim);

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
                contacts_.updateAccelerationDiff(data.contacts,
                                                 bottomRows(data.XAccx, get_ps().get_nv_dim()));
                contacts_.updateForceDiff(data.contacts, topRows(data.df_dx, nc_active_dim),
                                          topRows(data.df_du, nc_active_dim));
            }
            costs_.calcDiff(data.costs, x, u);
            if (constraints_.get_ng() > 0 || constraints_.get_nh() > 0)
            {
                constraints_.calcDiff(data.constraints, x, u);
            }
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            costs_.calcDiff(data.costs, x);
            if (constraints_.get_ng() > 0 || constraints_.get_nh() > 0)
            {
                constraints_.calcDiff(data.constraints, x);
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void quasiStatic(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlVectorType> &u,
                         const std::size_t maxiter, const NumScalar tol) const
        {
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, PS::DimNQ_t::Value> q =
                head(x, get_ps().get_nq_dim());
            auto nc_active_dim = contacts_.get_nc_active_dim();

            head(data.tmp_xstatic, get_ps().get_nq_dim()) = q;
            tail(data.tmp_xstatic, get_ps().get_nv_dim()).setZero();
            u.setZero();

            pinocchio::computeAllTerms(get_robot(), data.robot, q,
                                       tail(data.tmp_xstatic, get_ps().get_nv_dim()));
            pinocchio::computeJointJacobians(get_robot(), data.robot, q);
            pinocchio::rnea(get_robot(), data.robot, q,
                            tail(data.tmp_xstatic, get_ps().get_nv_dim()),
                            tail(data.tmp_xstatic, get_ps().get_nv_dim()));
            actuation_.calc(data.actuation, data.tmp_xstatic, u);
            actuation_.calcDiff(data.actuation, data.tmp_xstatic, u);
            contacts_.calc(data.contacts, data.tmp_xstatic);

            // Allocates memory
            data.tmp_Jstatic.conservativeResize(get_ps().get_nv(), get_ps().get_nu() + nc_active_dim.value());
            data.tmp_Jstatic.leftCols(get_ps().get_nu()) = data.actuation.dtau_du;
            rightCols(data.tmp_Jstatic, nc_active_dim) =
                topRows(data.contacts.Jc, nc_active_dim).transpose();
            u.noalias() = head(pseudoInverse(data.tmp_Jstatic) * data.robot.tau, get_ps().get_nu());
            data.robot.tau.setZero();
        }

        Data_t createData()
        {
            return Data_t(*this);
        }

        const CostModelManager_t &get_costs() const
        {
            return costs_;
        }

        const ConstraintModelManager_t &get_constraints() const
        {
            return constraints_;
        }

        const ContactModelManager_t &get_contacts() const
        {
            return contacts_;
        }

        const ActuationModel_t &get_actuation() const
        {
            return actuation_;
        }

        using Base::get_ps;

        using Base::get_robot;
        using Base::get_state;

    protected:
        CostModelManager_t costs_;
        ConstraintModelManager_t constraints_;
        ContactModelManager_t contacts_;
        ActuationModel_t actuation_;

        bool with_armature_;
        VectorNv_t armature_;
        NumScalar JMinvJt_damping_;
        bool enable_force_;

    }; // class NodeModelContactFwdDynTpl

} // namespace galileo

#endif // __galileo_predictive_nodes_node_contact_fwddyn_hpp__
