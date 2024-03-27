#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <mutex>

#include <pinocchio/fwd.hpp>

#include <Eigen/Core>

#include "galileo/legged-model/LeggedBody.h"
#include "galileo/legged-model/LeggedRobotProblemData.h"
#include "galileo/legged-model/LeggedRobotStates.h"
#include "galileo/legged-model/EnvironmentSurfaces.h"
#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/tools/GNUPlotInterface.h"
#include "galileo/tools/MeshcatInterface.h"

namespace galileo
{
    namespace legged
    {

        class LeggedInterface
        {
        public:
            using T_ROBOT_STATE = casadi::DM;
            using LeggedRobotProblemData = galileo::legged::constraints::LeggedRobotProblemData;
            using LeggedConstraintBuilderType =
                std::shared_ptr<galileo::opt::ConstraintBuilder<LeggedRobotProblemData>>;
            using EnvironmentSurfaces = galileo::legged::environment::EnvironmentSurfaces;
            using LeggedTrajOpt = galileo::opt::TrajectoryOpt<LeggedRobotProblemData, galileo::legged::contact::ContactMode>;

            LeggedInterface(std::string sol_data_dir = "../examples/visualization/solution_data/", std::string plot_dir = "../examples/visualization/plots/");

            /**
             * @brief Load the model from a file.
             */
            void LoadModel(std::string model_file_location, std::vector<std::string> end_effector_names);

            /**
             * @brief Load the parameters from a file.
             */
            void LoadParameters(std::string parameter_file_location);

            /**
             * @brief Create the contact sequence.
             */
            void setContactSequence(std::shared_ptr<contact::ContactSequence> contact_sequence)
            {
                assert(robot_ != nullptr);
                robot_->contact_sequence = contact_sequence;
                casadi::Dict empty_opts;
                robot_->fillModeDynamics(empty_opts);
                phases_set_ = true;
            }

            /**
             * @brief Create the contact sequence from basic data.
             */
            void setContactSequence(std::vector<int> knot_num, std::vector<double> knot_time, std::vector<uint> mask_vec,
                                    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces);

            /**
             * @brief Initialize the problem
             */
            void Initialize(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state);

            /**
             * @brief Solve the problem
             */
            void Update(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state);

            /**
             * @brief Get the solution.
             *
             * @param query_times The times at which to get the solution state and input.(num_times x 1 vector)
             * @param state_result The state at the query times (num_states x num_times)
             * @param input_result The input at the query times (num_inputs x num_times)
             */
            bool GetSolution(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result);

            /**
             * @brief Get the solution and plot the constraints
             */
            void VisualizeSolutionAndConstraints(const Eigen::VectorXd &query_times, Eigen::MatrixXd &state_result, Eigen::MatrixXd &input_result);

            /**
             * @brief Add a surface to the environment.
             */
            void addSurface(const environment::SurfaceData &surface) { surfaces_->push_back(surface); }

            /**
             * @brief Get the surfaces in the environment.
             */
            std::shared_ptr<EnvironmentSurfaces> &surfaces() { return surfaces_; }

            /**
             * @brief Get the legged constraint builders (velocity, friction cone, contact)
             */
            std::vector<LeggedConstraintBuilderType>
            getLeggedConstraintBuilders() const;

            /**
             * @brief Get the states of the robot.
             *
             * @return std::shared_ptr<LeggedRobotStates>
             */
            std::shared_ptr<LeggedRobotStates> states() const { return states_; }

            /**
             * @brief Get the trajectory optimizer.
             *
             * @return std::shared_ptr<LeggedTrajOpt>
             */
            std::shared_ptr<LeggedTrajOpt> getTrajectoryOptimizer() const { return trajectory_opt_; };

            /**
             * @brief Get the robot model.
             *
             * @return std::shared_ptr<LeggedBody>
             */
            std::shared_ptr<LeggedBody> getRobotModel() const { return robot_; };

            /**
             * @brief Get the end effectors of the robot.
             *
             */
            contact::RobotEndEffectors getEndEffectors() const { return robot_->getEndEffectors(); }

            /**
             * @brief Get the joint names in the order they are defined in the model.
             */
            std::vector<std::string> getJointNames() const { return robot_->model.names; }

            /**
             * @brief Get the horizon. This is not technically correct, it is actually the problem data dt.
             * The problem data can change after a solution is found, but it serves
             * as a best approximate for now.
             */
            double getSolutionDT()
            {
                assert(problem_data_ != nullptr);
                return problem_data_->phase_sequence->getDT();
            }

            // Checking solution state
            bool isRobotModelLoaded() const { return robot_ != nullptr; }
            bool isParametersLoaded() const { return parameters_set_; }
            bool isPhasesSet() const { return phases_set_; }
            bool isFullyInitialized() const { return fully_initialized_; }
            bool isSolutionSet() const { return solution_interface_ == nullptr ? false : solution_interface_->isSolutionSet(); }
            bool isSurfaceSet() const { return surfaces_ == nullptr ? false : surfaces_->size(); }
            bool CanInitialize() const { return isRobotModelLoaded() && isParametersLoaded() && isPhasesSet() && isSurfaceSet(); }

        private:
            /**
             * @brief Create the trajectory optimizer.
             */
            void CreateTrajOpt();

            /**
             * @brief Create the running and terminal costs.
             */
            void CreateCost(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state, casadi::Function &L, casadi::Function &Phi);

            /**
             * @brief update the problem data with new boundary conditions
             */
            void UpdateProblemBoundaries(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state);

            /**
             * @brief Create the problem data from the loaded parameter values
             */
            void CreateProblemData(const T_ROBOT_STATE &initial_state, const T_ROBOT_STATE &target_state);

            std::shared_ptr<LeggedRobotStates> states_; /**< Definition of the state. */

            std::shared_ptr<LeggedRobotProblemData> problem_data_; /**< The problem data. */

            std::shared_ptr<EnvironmentSurfaces> surfaces_; /**< The surfaces. */

            std::shared_ptr<contact::ContactSequence> contact_sequence_; /**< The gait. */

            std::shared_ptr<opt::solution::Solution> solution_interface_; /**< The solution interface. */
            std::mutex solution_mutex_;

            std::shared_ptr<tools::MeshcatInterface> meshcat_interface; /**< The meshcat interface. */

            std::shared_ptr<tools::GNUPlotInterface> plotting_interface; /**< The plotting interface. */

            casadi::Dict opts_;

            std::shared_ptr<opt::DecisionDataBuilder<LeggedRobotProblemData>> decision_builder_;

            std::shared_ptr<LeggedTrajOpt> trajectory_opt_; /**< The trajectory optimizer. */
            std::mutex trajectory_opt_mutex_;

            std::shared_ptr<LeggedBody> robot_; /**< The robot model. */

            std::string model_file_location_;

            struct CostParameters
            {
                Eigen::VectorXd Q_diag;
                Eigen::VectorXd R_diag;
            };
            CostParameters cost_params_;

            std::string body_name_ = "base";

            bool parameters_set_ = false;
            bool phases_set_ = false;

            bool fully_initialized_ = false;
        };

    }
}
