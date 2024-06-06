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

            std::string pub_orientation_rep_str;
            if (!nh_->getParam("/galileo_ros/published_orientation_representation", pub_orientation_rep_str))
            {
                /*See 'OrientationDefinition.h' for possible options*/
                ROS_INFO("No orientation representation provided, assuming quaternion [x, y, z, w]");
                published_orientation_representation_ = math::OrientationDefinition::Quaternion;
            }
            else
            {
                published_orientation_representation_ = static_cast<math::OrientationDefinition>(std::stoi(pub_orientation_rep_str));
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
            std::vector<contact::EE_Types> end_effector_types;
            for (size_t i = 0; i < msg->end_effector_names.size(); ++i)
            {
                end_effector_names.push_back(msg->end_effector_names[i]);
                end_effector_types.push_back(static_cast<contact::EE_Types>(msg->end_effector_types[i]));
            }

           internal_orientation_representation_ = static_cast<math::OrientationDefinition>(msg->orientation_definition);

            LoadModel(model_location, end_effector_names, end_effector_types, internal_orientation_representation_);
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

                contact_surface_ids.push_back(phase.mode.contact_surface_ids);
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

            Update(X0, Xf);
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

            Eigen::Index new_rows;
            if (published_orientation_representation_ == math::OrientationDefinition::Quaternion)
            {
                new_rows = states()->nh + states()->nqbp + 4 + states()->nvju;
                res.qj_index = states()->nh + states()->nqbp + 4;
            }
            else
            {
                new_rows = states()->nh + states()->nqbp + 3 + states()->nvju;
                res.qj_index = states()->nh + states()->nqbp + 3;
            }

            Eigen::MatrixXd new_state(new_rows, state_solution.cols());
            Eigen::MatrixXd new_orientation = galileo::math::OrientationConversion(state_solution.block(states()->qbo_index, 0, states()->nqbo, state_solution.cols()), internal_orientation_representation_, published_orientation_representation_);
            new_state << state_solution.topRows(states()->qbo_index), new_orientation, state_solution.bottomRows(states()->nvju);

            res.nx = new_rows;
            std::vector<double> state_vec = std::vector<double>(new_state.data(), new_state.data() + new_state.size());

            std::vector<double> input_vec(input_solution.data(), input_solution.data() + input_solution.size());

            res.joint_names = getJointNames();
            res.joint_names.erase(res.joint_names.begin(), res.joint_names.begin() + 2); // delete the universe and root_joint

            res.X_t_wrapped = state_vec;
            res.U_t_wrapped = input_vec;

            res.nu = input_solution.rows();

            res.solution_horizon = getSolutionDT();

            for (auto &ee : getEndEffectors())
            {
                res.end_effector_types.push_back(static_cast<int>(ee.second->ee_type));
            }

            auto contact_sequence = getRobotModel()->getContactSequence();

            galileo_ros::ContactSequence contact_sequence_msg;
            for (auto &phase : contact_sequence->getPhases())
            {
                galileo_ros::ContactPhase phase_msg;
                phase_msg.knot_num = phase.knot_points;
                phase_msg.knot_time = phase.time_value;

                phase_msg.mode.contact_surface_ids = phase.mode.contact_surfaces;
                phase_msg.mode.ee_names = getRobotModel()->getEndEffectorNames();

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