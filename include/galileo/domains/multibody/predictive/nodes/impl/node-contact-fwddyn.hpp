#ifndef __galileo_multibody_predictive_nodes_node_contact_fwddyn_hpp__
#define __galileo_multibody_predictive_nodes_node_contact_fwddyn_hpp__

#include <pinocchio/algorithm/centroidal.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/contact-dynamics.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/kinematics-derivatives.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <pinocchio/algorithm/rnea.hpp>

#include "galileo/common/math/matrix-decomposition.hpp"
#include "galileo/core/data/data-collector-default.hpp"
#include "galileo/domains/multibody/spatial/contacts/contact-manager.hpp"
#include "galileo/predictive/nodes/node-base.hpp"
#include "galileo/domains/multibody/predictive/nodes/fwd.hpp"

namespace galileo
{

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct NodeContactFwdDynTpl;

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct traits<NodeContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = NodeModelContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Data_t = NodeDataContactFwdDynTpl<PS, ContactCollectionTpl>;

        using DimNC_t = DimensionTpl<Eigen::Dynamic>;
        static constexpr int NC = DimNC_t::Value;
        using DimNU_t = typename PS::DimNUa_t;
        static constexpr int NU = DimNU_t::Value;

        using ContactManagerMeta_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using ContactModelManager_t = typename traits<ContactManagerMeta_t>::ModelManager_t;
        using ContactDataManager_t = typename traits<ContactManagerMeta_t>::DataManager_t;

        using DataCollector_t = DataCollectorContactTpl<PS, ContactCollectionTpl>;

        // using NC = traits<NodeDerived>::NC;
        // using Kinv_t = Eigen::GMatrix<VarScalar, PS::NV + NC, PS::NV + NC>;
        // using Jstatic_t = Eigen::GMatrix<VarScalar, PS::NV, NU + NC>;
        using Kinv_t =
            ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, AddDim_v<PS::NV, NC>, AddDim_v<PS::NV, NC>>>;
        using MatrixNcNdx_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NC, PS::NDX>>;
        using MatrixNcNu_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, NC, NU>>;
        using Jstatic_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, PS::NV, AddDim_v<NU, NC>>>;
    };

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct traits<NodeDataContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using Meta_t = NodeContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>;
    };

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct traits<NodeModelContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>>
    {
        using Meta_t = NodeContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>;
    };

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    struct NodeDataContactFwdDynTpl
        : public NodeDataBase<NodeDataContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = NodeDataBase<NodeDataContactFwdDynTpl<PS, ContactCollectionTpl>, PS>;

        using ContactManagerMeta_t = typename traits<Meta_t>::ContactManagerMeta_t;
        using ContactModelManager_t = typename traits<Meta_t>::ContactModelManager_t;
        using ContactDataManager_t = typename traits<Meta_t>::ContactDataManager_t;

        using XAcc_t = ArenaMatrixTpl<typename PS::XAcc_t>;
        using XAccx_t = ArenaMatrixTpl<typename PS::XAccx_t>;
        using XAccu_t = ArenaMatrixTpl<typename PS::XAccu_t>;
        using L_t = typename PS::L_t;
        using Lx_t = ArenaMatrixTpl<typename PS::Lx_t>;
        using Lu_t = ArenaMatrixTpl<typename PS::Lu_t>;
        using Lxx_t = ArenaMatrixTpl<typename PS::Lxx_t>;
        using Lxu_t = ArenaMatrixTpl<typename PS::Lxu_t>;
        using Luu_t = ArenaMatrixTpl<typename PS::Luu_t>;
        using H_t = ArenaMatrixTpl<typename PS::H_t>;
        using Hx_t = ArenaMatrixTpl<typename PS::Hx_t>;
        using Hu_t = ArenaMatrixTpl<typename PS::Hu_t>;
        using G_t = ArenaMatrixTpl<typename PS::G_t>;
        using Gx_t = ArenaMatrixTpl<typename PS::Gx_t>;
        using Gu_t = ArenaMatrixTpl<typename PS::Gu_t>;

        using CostDataManager_t = typename PS::CostDataManager_t;
        using ConstraintDataManager_t = typename PS::ConstraintDataManager_t;

        using DataCollector_t = typename traits<Meta_t>::DataCollector_t;
        using RobotData_t = typename DataCollector_t::RobotData_t;
        using ActuationData_t = typename DataCollector_t::ActuationData_t;
        using JointData_t = typename DataCollector_t::JointData_t;

        using VectorNx_t = ArenaMatrixTpl<typename PS::VectorNx_t>;

        using Kinv_t = typename traits<Meta_t>::Kinv_t;
        using MatrixNcNdx_t = typename traits<Meta_t>::MatrixNcNdx_t;
        using MatrixNcNu_t = typename traits<Meta_t>::MatrixNcNu_t;
        using Jstatic_t = typename traits<Meta_t>::Jstatic_t;

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
        G_t &G_accessor() { return constraints.H; }
        const G_t &G_accessor() const { return constraints.H; }
        Gx_t &Gx_accessor() { return constraints.Hx; }
        const Gx_t &Gx_accessor() const { return constraints.Hx; }
        Gu_t &Gu_accessor() { return constraints.Hu; }
        const Gu_t &Gu_accessor() const { return constraints.Hu; }

        NodeDataContactFwdDynTpl(const Model_t &model, MemoryArena &arena)
            : XAcc(arena, model.get_ps().get_nv()),
              XAccx(arena, model.get_ps().get_nv(), model.get_ps().get_ndx()),
              XAccu(arena, model.get_ps().get_nv(), model.get_ps().get_nu()),
              data_collector(std::make_shared<DataCollector_t>(
                  std::make_shared<RobotData_t>(model.get_state().get_robot()),
                  std::make_shared<ActuationData_t>(model.get_actuation().createData(arena)),
                  std::make_shared<JointData_t>(model.get_ps()))),
              robot(data_collector->robot),
              actuation(data_collector->actuation),
              joint(data_collector->joint),
              contacts(model.get_contacts().createData(arena, robot.get())),
              costs(model.get_costs().createData(arena, data_collector.get())),
              constraints(model.get_constraints().createData(arena, data_collector.get())),
              Kinv(arena,
                   model.get_ps().get_nv() + model.get_contacts().get_n_total(),
                   model.get_ps().get_nv() + model.get_contacts().get_n_total()),
              df_dx(arena, model.get_contacts().get_n_total(), model.get_ps().get_ndx()),
              df_du(arena, model.get_contacts().get_n_total(), model.get_ps().get_nu()),
              tmp_xstatic(arena, model.get_ps().get_nx()),
              tmp_Jstatic(arena, model.get_ps().get_nv(), model.get_ps().get_nu() + model.get_contacts().get_n_total())
        {
            XAcc.setZero();
            XAccx.setZero();
            XAccu.setZero();
            joint->dtau_du.diagonal().setOnes();
            Kinv.setZero();
            df_dx.setZero();
            df_du.setZero();
            tmp_xstatic.setZero();
            tmp_Jstatic.setZero();
            robot->lambda_c.resize(model.get_contacts().get_n_total());
            robot->lambda_c.setZero();
        }

        XAcc_t XAcc;
        XAccx_t XAccx;
        XAccu_t XAccu;

        std::shared_ptr<DataCollector_t> data_collector;
        std::shared_ptr<RobotData_t> robot;
        std::shared_ptr<ActuationData_t> actuation;
        std::shared_ptr<JointData_t> joint;
        ContactDataManager_t contacts;
        CostDataManager_t costs;
        ConstraintDataManager_t constraints;

        Kinv_t Kinv;
        MatrixNcNdx_t df_dx;
        MatrixNcNu_t df_du;
        VectorNx_t tmp_xstatic;
        Jstatic_t tmp_Jstatic;

    }; // class NodeDataContactFwdDynTpl

    template <typename PhaseSpec, template <typename> class ContactCollectionTpl>
    class NodeModelContactFwdDynTpl
        : public NodeModelBase<NodeModelContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = NodeModelBase<NodeModelContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>;

        using ContactManagerMeta_t = typename traits<Meta_t>::ContactManagerMeta_t;
        using ContactModelManager_t = typename traits<Meta_t>::ContactModelManager_t;
        using ContactDataManager_t = typename traits<Meta_t>::ContactDataManager_t;

        using CostModelManager_t = typename PS::CostModelManager_t;
        using ConstraintModelManager_t = typename PS::ConstraintModelManager_t;
        using State_t = typename PS::State_t;
        using ActuationModel_t = typename PS::ActuationModel_t;
        using NumScalar = typename PS::NumScalar;
        using VectorNv_t = typename PS::VectorNv_t;
        using MatrixX_t = typename PS::MatrixX_t;
        using VectorX_t = typename PS::VectorX_t;

        NodeModelContactFwdDynTpl(PS &ps,
                                  const State_t &state,
                                  const CostModelManager_t &costs,
                                  const ConstraintModelManager_t &constraints,
                                  const ContactModelManager_t &contacts,
                                  const ActuationModel_t &actuation,
                                  const NumScalar &JMinvJt_damping,
                                  const bool enable_force)
            : Base(ps, state),
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
            // get_ps().set_nu(nua);
            set_u_lb(NumScalar(-1.) * get_state().get_robot().effortLimit.tail(nua));
            set_u_ub(NumScalar(+1.) * get_state().get_robot().effortLimit.tail(nua));

            armature_.setZero();
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            const auto nc_dim = get_contacts().get_n_active_dim();

            const auto q = head(x, get_ps().get_nq_dim());
            const auto v = tail(x, get_ps().get_nv_dim());
            pinocchio::computeAllTerms(get_state().get_robot(), *data.robot, q, v);
            pinocchio::computeCentroidalMomentum(get_state().get_robot(), *data.robot);

            if (!with_armature_)
            {
                data.robot->M.diagonal() += armature_;
            }
            get_actuation().calc(*data.actuation.get(), x, u);
            get_contacts().calc(data.contacts, x);

            pinocchio::forwardDynamics(get_state().get_robot(),
                                       *data.robot,
                                       data.actuation->tau,
                                       topRows(data.contacts.Jc, nc_dim),
                                       head(data.contacts.a0, nc_dim),
                                       JMinvJt_damping_);
            data.XAcc = data.robot->ddq;

            get_contacts().updateAcceleration(data.contacts, data.robot->ddq);
            get_contacts().updateForce(data.contacts, data.robot->lambda_c);
            data.joint->a = data.robot->ddq;
            data.joint->tau = u;
            get_costs().calc(data.costs, x, u);
            data.L_accessor() = data.costs.L;
            if (get_constraints().get_n_active() > 0 || get_constraints().get_n_active() > 0)
            {
                // data.constraints.resize(this, data);
                get_constraints().calc(data.constraints, x, u);
            }
        }

        template <typename StateVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const auto q = head(x, get_ps().get_nq_dim());
            const auto v = tail(x, get_ps().get_nv_dim());
            pinocchio::computeAllTerms(get_state().get_robot(), *data.robot, q, v);
            pinocchio::computeCentroidalMomentum(get_state().get_robot(), *data.robot);
            get_costs().calc(data.costs, x);
            data.L_accessor() = data.costs.L;
            if (get_constraints().get_n_active() > 0 || get_constraints().get_n_active() > 0)
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
            const auto nv_dim = get_ps().get_nv_dim();
            const auto nc_dim = get_contacts().get_n_active_dim();
            const auto q = head(x, get_ps().get_nq_dim());
            const auto v = tail(x, get_ps().get_nv_dim());

            // Computing the dynamics derivatives
            // We resize the Kinv matrix because Eigen cannot call block operations
            // recursively: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=408. Therefore,
            // it is not possible to pass data.Kinv.topLeftCorner(nv + nc, nv + nc)
            data.Kinv.resize(get_ps().get_nv() + nc_dim.value(), nv_dim.value() + nc_dim.value());
            pinocchio::computeRNEADerivatives(
                get_state().get_robot(), *data.robot, q, v, data.XAcc, data.contacts.fext);
            get_contacts().updateRneaDiff(data.contacts, *data.robot);
            pinocchio::getKKTContactDynamicMatrixInverse(
                get_state().get_robot(), *data.robot, topRows(data.contacts.Jc, nc_dim), data.Kinv);

            get_actuation().calcDiff(*data.actuation.get(), x, u);
            get_contacts().calcDiff(data.contacts, x);

            const auto a_partial_dtau = topLeftCorner(data.Kinv, nv_dim, nv_dim);
            const auto a_partial_da = topRightCorner(data.Kinv, nv_dim, nc_dim);
            const auto f_partial_dtau = bottomLeftCorner(data.Kinv, nc_dim, nv_dim);
            const auto f_partial_da = bottomRightCorner(data.Kinv, nc_dim, nc_dim);

            leftCols(data.XAccx, nv_dim).noalias() = -a_partial_dtau * data.robot->dtau_dq;
            rightCols(data.XAccx, nv_dim).noalias() = -a_partial_dtau * data.robot->dtau_dv;
            data.XAccx.noalias() -= a_partial_da * topRows(data.contacts.da0_dx, nc_dim);
            data.XAccx.noalias() += a_partial_dtau * data.actuation->dtau_dx;
            data.XAccu.noalias() = a_partial_dtau * data.actuation->dtau_du;
            data.joint->da_dx = data.XAccx;
            data.joint->da_du = data.XAccu;

            // Computing the cost derivatives
            if (enable_force_)
            {
                topLeftCorner(data.df_dx, nc_dim, nv_dim).noalias() = f_partial_dtau * data.robot->dtau_dq;
                topRightCorner(data.df_dx, nc_dim, nv_dim).noalias() = f_partial_dtau * data.robot->dtau_dv;
                topRows(data.df_dx, nc_dim).noalias() += f_partial_da * topRows(data.contacts.da0_dx, nc_dim);
                topRows(data.df_dx, nc_dim).noalias() -= f_partial_dtau * data.actuation->dtau_dx;
                topRows(data.df_du, nc_dim).noalias() = -f_partial_dtau * data.actuation->dtau_du;
                get_contacts().updateAccelerationDiff(data.contacts, bottomRows(data.XAccx, nv_dim));
                get_contacts().updateForceDiff(data.contacts, topRows(data.df_dx, nc_dim), topRows(data.df_du, nc_dim));
            }
            get_costs().calcDiff(data.costs, x, u);
            if (get_constraints().get_n_active() > 0 || get_constraints().get_n_active() > 0)
            {
                get_constraints().calcDiff(data.constraints, x, u);
            }
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            get_costs().calcDiff(data.costs, x);
            if (get_constraints().get_n_active() > 0 || get_constraints().get_n_active() > 0)
            {
                get_constraints().calcDiff(data.constraints, x);
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void quasiStatic(Data_t &data,
                         const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlVectorType> &u,
                         const int maxiter,
                         const NumScalar tol) const
        {
            // Collect dimensions
            const auto nq_dim = get_ps().get_nq_dim();
            const auto nv_dim = get_ps().get_nv_dim();
            const auto nu_dim = get_ps().get_nu_dim();
            const auto nc_dim = get_contacts().get_n_active_dim();

            // Build the static state [ q; 0 ]
            const auto q = head(x, nq_dim);
            const auto v = VectorNv_t::Zero(nv_dim.value());
            head(data.tmp_xstatic, nq_dim) = q;
            tail(data.tmp_xstatic, nv_dim) = v;
            u.setZero();

            // Compute M(q) and bias h(q,0)
            pinocchio::computeAllTerms(get_state().get_robot(), *data.robot, q, v);
            pinocchio::rnea(get_state().get_robot(), *data.robot, q, v, v);
            const auto &h = data.robot->tau;

            // Linearize actuation & contacts
            get_actuation().calc(*data.actuation.get(), data.tmp_xstatic, u);
            get_actuation().calcDiff(*data.actuation.get(), data.tmp_xstatic, u);
            get_contacts().calc(data.contacts, data.tmp_xstatic);

            const auto B = data.actuation->dtau_du;
            const auto Jc = topRows(data.contacts.Jc, nc_dim);
            const auto a0 = head(data.contacts.a0, nc_dim);

            // data.tmp_Jstatic.conservativeResize(nv_dim.value(), nu_dim.value() + nc_dim.value());
            // leftCols(data.tmp_Jstatic, nu_dim) = B;
            // rightCols(data.tmp_Jstatic, nc_dim) = Jc.transpose();

            // // Solve [B  Jc^{T}] [u; \lambda] = h
            // VectorX_t z = pseudoInverse(data.tmp_Jstatic) * h;

            // data.robot->lambda_c = tail(z, nc_dim);
            // u = head(z, nu_dim);

            data.robot->tau.setZero();
        }

        Data_t createData(MemoryArena &arena) const { return Data_t(*this, arena); }

        const CostModelManager_t &get_costs() const { return costs_.get(); }
        const ConstraintModelManager_t &get_constraints() const { return constraints_.get(); }
        const ContactModelManager_t &get_contacts() const { return contacts_.get(); }
        const ActuationModel_t &get_actuation() const { return actuation_.get(); }

        using Base::get_ps;
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

#endif // __galileo_multibody_predictive_nodes_node_contact_fwddyn_hpp__
