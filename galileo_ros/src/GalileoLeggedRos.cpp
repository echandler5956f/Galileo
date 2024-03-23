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
            command_subscriber_ = nh_->subscribe(solver_id_ + "_command", 1, &GalileoLeggedRos::InitializationCallback, this);
        }

        void GalileoLeggedRos::InitServices()
        {
            init_state_service_ = nh_->advertiseService(solver_id_ + "_init_state_service", &GalileoLeggedRos::InitStateServiceCallback, this);
            get_solution_service_ = nh_->advertiseService(solver_id_ + "_get_solution", &GalileoLeggedRos::GetSolutionCallback, this);
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

                uint mask_bin = 0;
                for (auto &csurface_id : phase.contact_surface_ids)
                {

                    mask_bin << 1;
                    if (csurface_id != galileo::legged::environment::NO_SURFACE)
                    {
                        mask_bin |= 1;
                    }
                }
                mask_vec.push_back(mask_bin);
            }

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

            assert(X0.size() == states()->nx);
            assert(Xf.size() == states()->nx);

            Initialize(X0, Xf);
        }

        void GalileoLeggedRos::UpdateCallback(const galileo_ros::GalileoCommand::ConstPtr &msg)
        {
            auto X0 = casadi::DM(msg->X_initial);
            auto Xf = casadi::DM(msg->X_goal);

            assert(X0.size() == states()->nx);
            assert(Xf.size() == states()->nx);

            Update(X0, Xf);
        }

        bool GalileoLeggedRos::InitStateServiceCallback(galileo_ros::InitState::Request &req, galileo_ros::InitState::Response &res)
        {

            std::cout << "Checking if the solver can be initialized" << std::endl;
            res.model_set = RobotModelLoaded();
            res.solver_parameters_set = ParametersLoaded();
            res.contact_sequence_set = PhasesSet();
            res.fully_initted = FullyInitialized();
            res.solution_set = SolutionSet();

            return true;
        }

        bool GalileoLeggedRos::GetSolutionCallback(galileo_ros::SolutionRequest::Request &req, galileo_ros::SolutionRequest::Response &res)
        {
            Eigen::MatrixXd state_solution;
            Eigen::MatrixXd input_solution;

            Eigen::Map<Eigen::VectorXd> query_times(req.times.data(), req.times.size());
            bool solution_exists = GetSolution(query_times, state_solution, input_solution);

            if (!solution_exists)
            {
                res.solution_exists = false;
                return false;
            }

            std::vector<double> state_vec(state_solution.data(), state_solution.data() + state_solution.size());
            std::vector<double> input_vec(input_solution.data(), input_solution.data() + input_solution.size());

            res.joint_names = getJointNames();
            res.qj_index = states()->qj_index;

            res.X_t_wrapped = state_vec;
            res.U_t_wrapped = input_vec;

            res.nx = state_solution.rows();
            res.nu = input_solution.rows();

            res.solution_horizon = getSolutionDT();

            res.times_evaluated = req.times;
            res.solution_exists = true;

            return true;
        }
    }
}