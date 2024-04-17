#include "galileo_ros/GalileoLeggedRos.h"

namespace galileo
{
    namespace legged
    {
        void GalileoLeggedRos::InitPublishers() {}

        void GalileoLeggedRos::InitSubscribers()
        {
            // Subscribe to the model location
            model_location_subscriber_ = nh_->subscribe(solver_id_ + "_model_location", 1, &GalileoLeggedRos::ModelLocationCallback, this);

            // Subscribe to the parameter location
            parameter_location_subscriber_ = nh_->subscribe(solver_id_ + "_parameter_location", 1, &GalileoLeggedRos::ParameterLocationCallback, this);

            // Subscribe to the contact sequence
            contact_sequence_subscriber_ = nh_->subscribe(solver_id_ + "_contact_sequence", 1000, &GalileoLeggedRos::ContactSequenceCallback, this);

            // Subscribe to the environment surface
            surface_subscriber_ = nh_->subscribe(solver_id_ + "_add_environment_surface", 1000, &GalileoLeggedRos::SurfaceCallback, this);

            // Subscribe to the initialization command
            command_subscriber_ = nh_->subscribe(solver_id_ + "_command", 1, &GalileoLeggedRos::GeneralCommandCallback, this);

            if (!nh_->getParam("galileo_ros/orientation_definition", orientation_definition_))
            {
                /*  This parameter ensures that the orientation from the solution is converted to the desired format. By default (-1), the orientation is described as a quaternion [x, y, z, w]

                    The following are the possible orientation definitions:
                        zyx = 0
                        zyz = 1
                        zxy = 2
                        zxz = 3
                        yxz = 4
                        yxy = 5
                        yzx = 6
                        yzy = 7
                        xyz = 8
                        xyx = 9
                        xzy = 10
                        xzx = 11
                    From galileo/math/Quat2Euler.h
                */
                ROS_INFO("No orientation definition provided, assuming quaternion [x, y, z, w]");
            }
        }

        void GalileoLeggedRos::InitServices()
        {
            init_state_service_ = nh_->advertiseService(solver_id_ + "_init_state_service", &GalileoLeggedRos::InitStateServiceCallback, this);
            get_solution_service_ = nh_->advertiseService(solver_id_ + "_get_solution", &GalileoLeggedRos::GetSolutionCallback, this);
            write_sol_plot_cons_service_ = nh_->advertiseService(solver_id_ + "_write_sol_plot_cons", &GalileoLeggedRos::WriteSolPlotConsCallback, this);
        }

        void GalileoLeggedRos::ModelLocationCallback(const galileo_ros::RobotModel::ConstPtr &msg)
        {
            std::string model_location = msg->model_file_location;
            std::vector<std::string> end_effector_names;
            for (size_t i = 0; i < msg->end_effector_names.size(); ++i)
            {
                end_effector_names.push_back(msg->end_effector_names[i]);
            }

            LoadModel(model_location, end_effector_names);
        }

        void GalileoLeggedRos::ParameterLocationCallback(const galileo_ros::ParameterFileLocation::ConstPtr &msg)
        {
            std::string parameter_location = msg->parameter_file_location;
            LoadParameters(parameter_location);
        }

        void GalileoLeggedRos::ContactSequenceCallback(const galileo_ros::ContactSequence::ConstPtr &msg)
        {
            std::vector<uint> mask_vec;
            std::vector<int> knot_nums;
            std::vector<double> knot_times;

            std::vector<std::vector<environment::SurfaceID>> contact_surface_ids;
            for (auto &phase : msg->phases)
            {
                knot_nums.push_back(phase.knot_num);
                knot_times.push_back(phase.knot_time);

                contact_surface_ids.push_back(phase.contact_surface_ids);
            }

            mask_vec = helper::getMaskVectorFromContactSurfaces(contact_surface_ids);

            setContactSequence(knot_nums, knot_times, mask_vec, contact_surface_ids);
        }

        void GalileoLeggedRos::SurfaceCallback(const galileo_ros::EnvironmentSurface::ConstPtr &msg)
        {
            // UNIMPLEMENTED
            addSurface(environment::createInfiniteGround());
        }

        void GalileoLeggedRos::InitializationCallback(const galileo_ros::GalileoCommand::ConstPtr &msg)
        {
            auto X0 = casadi::DM(msg->X_initial);
            auto Xf = casadi::DM(msg->X_goal);

            assert(X0.size1() == states()->nx);
            assert(Xf.size1() == states()->nx);

            Initialize(X0, Xf);
        }

        void GalileoLeggedRos::UpdateCallback(const galileo_ros::GalileoCommand::ConstPtr &msg)
        {
            auto X0 = casadi::DM(msg->X_initial);
            auto Xf = casadi::DM(msg->X_goal);

            assert(X0.size1() == states()->nx);
            assert(Xf.size1() == states()->nx);

            Update(0., X0, Xf);
        }

        bool GalileoLeggedRos::InitStateServiceCallback(galileo_ros::InitState::Request &req, galileo_ros::InitState::Response &res)
        {

            res.model_set = isRobotModelLoaded();

            res.solver_parameters_set = isParametersLoaded();

            res.contact_sequence_set = isPhasesSet();

            res.environment_surface_set = isSurfaceSet();

            res.fully_initted = isFullyInitialized();

            res.solution_set = isSolutionSet();

            return true;
        }

        void GalileoLeggedRos::GeneralCommandCallback(const galileo_ros::GalileoCommand::ConstPtr &msg)
        {
            if (msg->command_type == "init")
            {
                InitializationCallback(msg);
            }

            UpdateCallback(msg);
        }

        bool GalileoLeggedRos::GetSolutionCallback(galileo_ros::SolutionRequest::Request &req, galileo_ros::SolutionRequest::Response &res)
        {
            Eigen::MatrixXd state_solution = Eigen::MatrixXd::Zero(states()->nx, req.times.size());
            Eigen::MatrixXd input_solution = Eigen::MatrixXd::Zero(states()->nu, req.times.size());

            Eigen::Map<Eigen::VectorXd> query_times(req.times.data(), req.times.size());

            bool solution_exists = GetSolution(query_times, state_solution, input_solution);

            if (!solution_exists)
            {
                res.solution_exists = false;
                return false;
            }
            std::vector<double> state_vec;
            if (orientation_definition_ != -1)
            {
                Eigen::MatrixXd new_state(state_solution.rows() - 1, state_solution.cols());
                Eigen::MatrixXd euler_angles = galileo::math::quaternion2Euler(state_solution.block(states()->q_index + 3, 0, 4, state_solution.cols()), orientation_definition_);
                new_state << state_solution.topRows(states()->q_index + 3), euler_angles, state_solution.bottomRows(states()->nq - states()->nqb);

                res.nx = state_solution.rows() - 1;
                state_vec = std::vector<double>(new_state.data(), new_state.data() + new_state.size());
            }
            else
            {
                res.nx = state_solution.rows();
                state_vec = std::vector<double>(state_solution.data(), state_solution.data() + state_solution.size());
            }

            std::vector<double> input_vec(input_solution.data(), input_solution.data() + input_solution.size());

            res.joint_names = getJointNames();
            res.joint_names.erase(res.joint_names.begin(), res.joint_names.begin() + 2); // delete the universe and root_joint

            res.qj_index = states()->qj_index;

            res.X_t_wrapped = state_vec;
            res.U_t_wrapped = input_vec;

            res.nu = input_solution.rows();

            res.solution_horizon = getSolutionDT();

            std::vector<int> ee_dofs;

            for (auto &ee : getEndEffectors())
            {
                ee_dofs.push_back(ee.second->is_6d ? 6 : 3);
            }

            res.ee_dofs = ee_dofs;

            // TODO: Add the contact modes to the service response so the WBC knows which end effectors are in contact

            auto contact_sequence = robot_->getContactSequence();

            galileo_ros::ContactSequence contact_sequence_msg;
            for (auto &phase : contact_sequence->getPhases())
            {
                galileo_ros::ContactPhase phase_msg;
                phase_msg.knot_num = phase.knot_points;
                phase_msg.knot_time = phase.time_value;

                for (auto &contact_surface_id : phase.mode.contact_surfaces)
                {
                    phase_msg.contact_surface_ids.push_back(contact_surface_id);
                }

                contact_sequence_msg.phases.push_back(phase_msg);
            }

            res.times_evaluated = req.times;
            res.solution_exists = true;

            return true;
        }

        bool GalileoLeggedRos::WriteSolPlotConsCallback(galileo_ros::WriteSolPlotCons::Request &req, galileo_ros::WriteSolPlotCons::Response &res)
        {
            Eigen::MatrixXd state_solution = Eigen::MatrixXd::Zero(states()->nx, req.times.size());
            Eigen::MatrixXd input_solution = Eigen::MatrixXd::Zero(states()->nu, req.times.size());

            Eigen::Map<Eigen::VectorXd> query_times(req.times.data(), req.times.size());
            VisualizeSolutionAndConstraints(query_times, state_solution, input_solution);

            return true;
        }
    }
}