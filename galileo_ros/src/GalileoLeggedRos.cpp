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
            can_init_service_ = nh_->advertiseService(solver_id_ + "_can_init_service", &GalileoLeggedRos::CanInitServiceCallback, this);
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

            casadi::MXVector solution;
            Update(X0, Xf, solution);
        }

        bool GalileoLeggedRos::CanInitServiceCallback(std_srvs::Trigger::Request &req, std_srvs::Trigger::Response &res)
        {
            res.success = CanInitialize();
            return true;
        }
    }
}