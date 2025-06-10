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

#include "galileo/multibody/contacts/fwd.hpp"
#include "galileo/multibody/contacts/contact-manager.hpp"

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

        // TODO: Refactor to also be able to handle when NC is known at compile time
        // static constexpr int NC = ContactDataManager_t::NC;
        static constexpr int NC = Eigen::Dynamic;
        static constexpr int NU = PS::NUa;

        using ContactManagerMeta_t = ContactManagerTpl<PS, ContactCollectionTpl>;
        using ContactModelManager_t = typename traits<ContactManagerMeta_t>::ModelManager_t;
        using ContactDataManager_t = typename traits<ContactManagerMeta_t>::DataManager_t;

        // using NC = traits<NodeDerived>::NC;
        // using Kinv_t = Eigen::Matrix<VarScalar, PS::NV + NC, PS::NV + NC>;
        // using Jstatic_t = Eigen::Matrix<VarScalar, PS::NV, NU + NC>;
        using Kinv_t = Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, Eigen::Dynamic>;
        using MatrixNcNdx_t = Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, PS::NDX>;
        using MatrixNcNu_t = Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, NU>;
        using Jstatic_t = Eigen::Matrix<typename PS::VarScalar, PS::NV, Eigen::Dynamic>;

        using MatrixNvNc_t = Eigen::Matrix<typename PS::VarScalar, PS::NV, Eigen::Dynamic>;
        using MatrixNcNv_t = Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, PS::NV>;
        using MatrixNc_t = Eigen::Matrix<typename PS::VarScalar, Eigen::Dynamic, Eigen::Dynamic>;
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
    struct NodeDataContactFwdDynTpl : public NodeDataBase<NodeDataContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>
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

        using Kinv_t = typename traits<Meta_t>::Kinv_t;
        using MatrixNcNdx_t = typename traits<Meta_t>::MatrixNcNdx_t;
        using MatrixNcNu_t = typename traits<Meta_t>::MatrixNcNu_t;
        using Jstatic_t = typename traits<Meta_t>::Jstatic_t;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        DEFAULT_ACCESSOR(ActuationData_t, actuation);
        DEFAULT_ACCESSOR(ConstraintDataManager_t, constraints);
        DEFAULT_ACCESSOR(CostDataManager_t, costs);

        DEFAULT_ACCESSOR(XAcc_t, XAcc);
        DEFAULT_ACCESSOR(XAccx_t, XAccx);
        DEFAULT_ACCESSOR(XAccu_t, XAccu);

        DEFAULT_ACCESSOR(L_t, L);
        DEFAULT_ACCESSOR(Lx_t, Lx);
        DEFAULT_ACCESSOR(Lu_t, Lu);
        DEFAULT_ACCESSOR(Lxx_t, Lxx);
        DEFAULT_ACCESSOR(Lxu_t, Lxu);
        DEFAULT_ACCESSOR(Luu_t, Luu);

        DEFAULT_ACCESSOR(H_t, H);
        DEFAULT_ACCESSOR(Hx_t, Hx);
        DEFAULT_ACCESSOR(Hu_t, Hu);

        DEFAULT_ACCESSOR(G_t, G);
        DEFAULT_ACCESSOR(Gx_t, Gx);
        DEFAULT_ACCESSOR(Gu_t, Gu);

        ContactDataManager_t contacts;
        Kinv_t Kinv;
        MatrixNcNdx_t df_dx;
        MatrixNcNu_t df_du;
        VectorNx_t tmp_xstatic;
        Jstatic_t tmp_Jstatic;

        ActuationData_t actuation;
        ConstraintDataManager_t constraints;
        CostDataManager_t costs;
        RobotData_t robot;

        XAcc_t XAcc;
        XAccx_t XAccx;
        XAccu_t XAccu;

        L_t L;
        Lx_t Lx;
        Lu_t Lu;
        Lxx_t Lxx;
        Lxu_t Lxu;
        Luu_t Luu;

        H_t H;
        Hx_t Hx;
        Hu_t Hu;

        G_t G;
        Gx_t Gx;
        Gu_t Gu;

        // template <typename DataCollector>
        // NodeDataContactFwdDynTpl(const Model_t &model, DataCollector *const collector)
        //     : XAcc(NV),
        //       XAccx(NV, NDX),
        //       XAccu(NV, model->get_nu()),
        //       L(VarScalar(0.)),
        //       r(model->get_nr()),
        //       Lx(NDX),
        //       Lu(model->get_nu()),
        //       Lxx(NDX, NDX),
        //       Lxu(NDX, model->get_nu()),
        //       Luu(model->get_nu(), model->get_nu()),
        //       H(model->get_nh()),
        //       Hx(model->get_nh(), NDX),
        //       Hu(model->get_nh(), model->get_nu()),
        //       G(model->get_ng()),
        //       Gx(model->get_ng(), NDX),
        //       Gu(model->get_ng(), model->get_nu()),
        //       pinocchio(pinocchio::DataTpl<Scalar>(model->get_pinocchio())),
        //       multibody(
        //           &pinocchio, model->get_actuation()->createData(),
        //           boost::make_shared<JointDataAbstract>(
        //               model->get_state(), model->get_actuation(), model->get_nu()),
        //           model->get_contacts()->createData(&pinocchio)),
        //       costs(model->get_costs()->createData(&multibody)),
        //       Kinv(NV +
        //                model->get_contacts()->get_nc_total(),
        //            NV +
        //                model->get_contacts()->get_nc_total()),
        //       df_dx(model->get_contacts()->get_nc_total(),
        //             NDX),
        //       df_du(model->get_contacts()->get_nc_total(), model->get_nu()),
        //       tmp_xstatic(model->get_state()->get_nx()),
        //       tmp_Jstatic(NV,
        //                   model->get_nu() + model->get_contacts()->get_nc_total())
        // {
        //     XAcc.setZero();
        //     XAccx.setZero();
        //     XAccu.setZero();
        //     r.setZero();
        //     Lx.setZero();
        //     Lu.setZero();
        //     Lxx.setZero();
        //     Lxu.setZero();
        //     Luu.setZero();
        //     H.setZero();
        //     Hx.setZero();
        //     Hu.setZero();
        //     G.setZero();
        //     Gx.setZero();
        //     Gu.setZero();
        //     multibody.joint->dtau_du.diagonal().setOnes();
        //     costs->shareMemory(this);
        //     if (model->get_constraints() != nullptr)
        //     {
        //         constraints = model->get_constraints()->createData(&multibody);
        //         constraints->shareMemory(this);
        //     }
        //     Kinv.setZero();
        //     df_dx.setZero();
        //     df_du.setZero();
        //     tmp_xstatic.setZero();
        //     tmp_Jstatic.setZero();
        //     pinocchio.lambda_c.resize(model->get_contacts()->get_nc_total());
        //     pinocchio.lambda_c.setZero();
        // }

    }; // class NodeDataFreeFwdTpl

    template <typename PhaseSpec,
              template <typename PS> class ContactCollectionTpl>
    class NodeModelContactFwdDynTpl : public NodeModelBase<NodeModelContactFwdDynTpl<PhaseSpec, ContactCollectionTpl>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = NodeContactFwdDynTpl<PS, ContactCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

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
            const std::size_t nc = contacts_->nc();

            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, NQ> q =
                x.head(PS::NQ);
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, NV> v =
                x.tail(PS::NV);

            // Computing the forward dynamics with the holonomic constraints defined by
            // the contact model
            pinocchio::computeAllTerms(*robot_, data.robot, q, v);
            pinocchio::computeCentroidalMomentum(*robot_, data.robot);

            if (!with_armature_)
            {
                data.robot.M.diagonal() += armature_;
            }
            actuation_->calc(data.multibody.actuation, x, u);
            contacts_->calc(data.multibody.contacts, x);

            pinocchio::forwardDynamics(
                *robot_, data.robot, data.multibody.actuation.tau,
                data.multibody.contacts.Jc.topRows(nc), data.multibody.contacts.a0.head(nc),
                JMinvJt_damping_);
            data.XAcc = data.robot.ddq;
            contacts_->updateAcceleration(data.multibody.contacts, data.robot.ddq);
            contacts_->updateForce(data.multibody.contacts, data.robot.lambda_c);
            data.multibody.joint.a = data.robot.ddq;
            data.multibody.joint.tau = u;
            costs_->calc(data.costs, x, u);
            data.cost = data.costs.cost;
            if (constraints_->ng() > 0 || constraints_->nh() > 0)
            {
                data.constraints.resize(this, data);
                constraints_->calc(data.constraints, x, u);
            }
        }

        template <typename StateVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x) const
        {
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, NQ> q =
                x.head(PS::NQ);
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, NV> v =
                x.tail(PS::NV);

            pinocchio::computeAllTerms(*robot_, data.robot, q, v);
            pinocchio::computeCentroidalMomentum(*robot_, data.robot);
            costs_->calc(data.costs, x);
            data.cost = data.costs.cost;
            if (constraints_->ng() > 0 || constraints_->nh() > 0)
            {
                data.constraints.resize(this, data);
                constraints_->calc(data.constraints, x);
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            const std::size_t nc = contacts_->nc();
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, NQ> q =
                x.head(PS::NQ);
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, NV> v =
                x.tail(PS::NV);

            // Computing the dynamics derivatives
            // We resize the Kinv matrix because Eigen cannot call block operations
            // recursively: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=408. Therefore,
            // it is not possible to pass data.Kinv.topLeftCorner(nv + nc, nv + nc)
            data.Kinv.resize(PS::NV + nc, PS::NV + nc);
            pinocchio::computeRNEADerivatives(*robot_, data.robot, q, v, data.XAcc,
                                              data.multibody.contacts.fext);
            contacts_->updateRneaDiff(data.multibody.contacts, data.robot);
            pinocchio::getKKTContactDynamicMatrixInverse(
                *robot_, data.robot, data.multibody.contacts.Jc.topRows(nc), data.Kinv);

            actuation_->calcDiff(data.multibody.actuation, x, u);
            contacts_->calcDiff(data.multibody.contacts, x);

            const Eigen::Block<MatrixNv_t> a_partial_dtau = data.Kinv.topLeftCorner(PS::NV, PS::NV);
            const Eigen::Block<MatrixNvNc_t> a_partial_da = data.Kinv.topRightCorner(PS::NV, nc);
            const Eigen::Block<MatrixNcNv_t> f_partial_dtau = data.Kinv.bottomLeftCorner(nc, PS::NV);
            const Eigen::Block<MatrixNc_t> f_partial_da = data.Kinv.bottomRightCorner(nc, nc);

            data.XAccx.leftCols(PS::NV).noalias() = -a_partial_dtau * data.robot.dtau_dq;
            data.XAccx.rightCols(PS::NV).noalias() = -a_partial_dtau * data.robot.dtau_dv;
            data.XAccx.noalias() -= a_partial_da * data.multibody.contacts.da0_dx.topRows(nc);
            data.XAccx.noalias() += a_partial_dtau * data.multibody.actuation.dtau_dx;
            data.XAccu.noalias() = a_partial_dtau * data.multibody.actuation.dtau_du;
            data.multibody.joint.da_dx = data.XAccx;
            data.multibody.joint.da_du = data.XAccu;

            // Computing the cost derivatives
            if (enable_force_)
            {
                data.df_dx.topLeftCorner(nc, PS::NV).noalias() =
                    f_partial_dtau * data.robot.dtau_dq;
                data.df_dx.topRightCorner(nc, PS::NV).noalias() =
                    f_partial_dtau * data.robot.dtau_dv;
                data.df_dx.topRows(nc).noalias() +=
                    f_partial_da * data.multibody.contacts->da0_dx.topRows(nc);
                data.df_dx.topRows(nc).noalias() -=
                    f_partial_dtau * data.multibody.actuation->dtau_dx;
                data.df_du.topRows(nc).noalias() =
                    -f_partial_dtau * data.multibody.actuation->dtau_du;
                contacts_->updateAccelerationDiff(data.multibody.contacts,
                                                  data.XAccx.bottomRows(PS::NV));
                contacts_->updateForceDiff(data.multibody.contacts, data.df_dx.topRows(nc),
                                           data.df_du.topRows(nc));
            }
            costs_->calcDiff(data.costs, x, u);
            if (constraints_->ng() > 0 || constraints_->nh() > 0)
            {
                constraints_->calcDiff(data.constraints, x, u);
            }
        }

        template <typename StateVectorType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x) const
        {
            costs_->calcDiff(data.costs, x);
            if (constraints_->ng() > 0 || constraints_->nh() > 0)
            {
                constraints_->calcDiff(data.constraints, x);
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void quasiStatic(Data_t &data, const Eigen::MatrixBase<StateVectorType> &x,
                         Eigen::MatrixBase<ControlVectorType> &u,
                         const std::size_t maxiter, const NumScalar tol) const
        {
            const Eigen::VectorBlock<const Eigen::Ref<const VectorNx_t>, NQ> q =
                x.head(PS::NQ);
            const std::size_t nc = contacts_.nc();

            data.tmp_xstatic.head(PS::NQ) = q;
            data.tmp_xstatic.tail(PS::NV).setZero();
            u.setZero();

            pinocchio::computeAllTerms(*robot_, data.robot, q,
                                       data.tmp_xstatic.tail(PS::NV));
            pinocchio::computeJointJacobians(*robot_, data.robot, q);
            pinocchio::rnea(*robot_, data.robot, q, data.tmp_xstatic.tail(PS::NV),
                            data.tmp_xstatic.tail(PS::NV));
            actuation_->calc(data.multibody.actuation, data.tmp_xstatic, u);
            actuation_->calcDiff(data.multibody.actuation, data.tmp_xstatic, u);
            contacts_->calc(data.multibody.contacts, data.tmp_xstatic);

            // Allocates memory
            data.tmp_Jstatic.conservativeResize(PS::NV, PS::NU + nc);
            data.tmp_Jstatic.leftCols(PS::NU) = data.multibody.actuation.dtau_du;
            data.tmp_Jstatic.rightCols(nc) =
                data.multibody.contacts.Jc.topRows(nc).transpose();
            u.noalias() = (pseudoInverse(data.tmp_Jstatic) * data.robot.tau).head(PS::NU);
            data.robot.tau.setZero();
        }

    protected:
        State_t *state_;
        ActuationModel_t *actuation_;
        CostModelManager_t *costs_;
        ConstraintModelManager_t *constraints_;
        ContactModelManager_t *contacts_;

        RobotModel_t *robot_;
        bool with_armature_;
        VectorNv_t armature_;
        NumScalar JMinvJt_damping_;
        bool enable_force_;

    }; // class NodeModelContactFwdDynTpl

} // namespace galileo

#endif // __galileo_predictive_nodes_node_contact_fwddyn_hpp__