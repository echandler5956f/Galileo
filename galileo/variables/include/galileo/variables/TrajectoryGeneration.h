#pragma once

#include "galileo/variables/States.h"
#include "galileo/model/ContactSequence.h"

namespace galileo
{
    namespace variables
    {
        /**
         * @brief 
         * 
         * @tparam Sym 
         */
        template <class Sym = Eigen::VectorXd>
        struct Target
        {
            /**
             * @brief Construct a new Target object
             * 
             * @param target_vars_ 
             * @param state_def_ 
             */
            Target(Sym target_vars_, States* state_def_)
            {
                this->target_vars = target_vars_;
                this->state_def = state_def_;
                Q = Eigen::MatrixXd::Identity(this->state_def->ndx, this->state_def->ndx);
            }

            /**
             * @brief 
             * 
             */
            Sym target_vars;

            /**
             * @brief 
             * 
             */
            States* state_def;
            
            /**
             * @brief Cost matrix
             * 
             */
            Eigen::MatrixXd Q;
        };

        /**
         * @brief 
         * 
         * @tparam Sym 
         */
        template <class Sym = Eigen::VectorXd>
        struct InitialCondition
        {
            /**
             * @brief Construct a new Initial Condition object
             * 
             * @param x0_vars_ 
             * @param state_def_ 
             * @param initial_contacts_ 
             */
            InitialCondition(Sym x0_vars_, States* state_def_, contact::ContactMode initial_contacts_)
            {
                this->x0_vars = x0_vars_;
                this->state_def = state_def_;
                this->initial_contacts = initial_contacts_;
            }

            /**
             * @brief 
             * 
             */
            Sym x0_vars;

            /**
             * @brief 
             * 
             */
            States* state_def;

            /**
             * @brief 
             * 
             */
            contact::ContactMode initial_contacts;
        };

        /**
         * @brief 
         * 
         * @tparam Sym 
         */
        template <class Sym = Eigen::VectorXd>
        struct ProblemSetup
        {
            /**
             * @brief Construct a new Problem Setup object
             * 
             * @param initial_condition_ 
             * @param contact_sequence_ 
             */
            ProblemSetup(InitialCondition<Sym> initial_condition_, contact::ContactSequence *contact_sequence_) : 
                initial_condition(initial_condition_),
                contact_sequence(contact_sequence_){}

            /**
             * @brief 
             * 
             * @return true 
             * @return false 
             */
            bool CheckValidity();

            /**
             * @brief 
             * 
             */
            InitialCondition<Sym> initial_condition;

            /**
             * @brief 
             * 
             */
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