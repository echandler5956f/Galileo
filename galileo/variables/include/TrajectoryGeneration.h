#include "States.h"
#include "../../model/include/ContactSequence.h"

namespace acro
{
    namespace variables
    {
        template <class Sym = Eigen::VectorXd>
        struct Target
        {
            Target(Sym target_vars_, States* state_def_)
            {
                this->target_vars = target_vars_;
                this->state_def = state_def_;
                Q = Eigen::MatrixXd::Identity(this->state_def->ndx, this->state_def->ndx);
            }

            Sym target_vars;
            States* state_def;
            // Cost matrix
            Eigen::MatrixXd Q;
        };

        template <class Sym = Eigen::VectorXd>
        struct InitialCondition
        {
            InitialCondition(Sym x0_vars_, States* state_def_, contact::ContactMode initial_contacts_)
            {
                this->x0_vars = x0_vars_;
                this->state_def = state_def_;
                this->initial_contacts = initial_contacts_;
            }

            Sym x0_vars;
            States* state_def;
            contact::ContactMode initial_contacts;
        };

        template <class Sym = Eigen::VectorXd>
        struct ProblemSetup
        {
            ProblemSetup(InitialCondition<Sym> initial_condition_, contact::ContactSequence *contact_sequence_) : 
                initial_condition(initial_condition_),
                contact_sequence(contact_sequence_){}

            bool CheckValidity();

            InitialCondition<Sym> initial_condition;
            contact::ContactSequence *contact_sequence;
        };

        // struct ProblemData
        // {
        //     std::shared_ptr<pinocchio::Model> model;
        //     std::shared_ptr<environment::EnvironmentSurfaces> environment_surfaces;
        //     contact::ContactSequence contact_sequence;
        //     // footstep trajectory parameters or evaluator
        //     std::vector<int> collocation_points_per_knot;
        //     std::vector<std::shared_ptr<variables::ConstraintBuilder>> constraint_builders;
        // };
    }
}