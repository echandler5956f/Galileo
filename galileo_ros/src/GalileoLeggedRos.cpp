#include "galileo_ros/GalileoLeggedRos.h"

namespace galileo
{
    namespace legged
    {
        void GalileoLeggedRos::InitPublishers() {}

        void GalileoLeggedRos::InitSubscribers()
        {
            // Subscribe to the model location
            model_location_subscriber_ = nh_->subscribe(solver_id + "_model_location", 1, &GalileoLeggedRos::ModelLocationCallback, this);

            // Subscribe to the parameter location
            parameter_location_subscriber_ = nh_->subscribe(solver_id + "_parameter_location", 1, &GalileoLeggedRos::ParameterLocationCallback, this);

            // Subscribe to the contact sequence
            contact_sequence_subscriber_ = nh_->subscribe(solver_id + "_contact_sequence", 1000, &GalileoLeggedRos::ContactSequenceCallback, this);

            // Subscribe to the environment surface
            surface_subscriber_ = nh_->subscribe(solver_id + "_add_environment_surface", 1000, &GalileoLeggedRos::SurfaceCallback, this);

            // Subscribe to the initialization command
            command_subscriber_ = nh_->subscribe(solver_id + "_command", 1, &GalileoLeggedRos::InitializationCallback, this);
        }

        void GalileoLeggedRos::InitServices()
        {
            can_init_service_ = nh_->advertiseService(solver_id + "_can_init_service", &GalileoLeggedRos::CanInitServiceCallback, this);
        }

        void GalileoLeggedRos::ModelLocationCallback(const galileo_ros::msg::ModelLocation::ConstPtr &msg)
        {
            std::string model_location = msg->model_location;
            std::vector<std::string> end_effector_names;
            for (size_t i = 0; i < msg->end_effector_names.size(); ++i)
            {
                end_effector_names.push_back(msg->end_effector_names[i]);
            }

            LoadModel(model_location, end_effector_names);
        }

        void GalileoLeggedRos::ParameterLocationCallback(const galileo_ros::msg::ParameterLocation::ConstPtr &msg)
        {
            std::string parameter_location = msg->parameter_location;
            LoadParameters(parameter_location);
        }

        void GalileoLeggedRos::ContactSequenceCallback(const galileo_ros::msg::ContactSequence::ConstPtr &msg)
        {
            std::vector<uint> mask_vec;

            int num_knot_nums = msg->knot_nums.size();
            assert(num_knot_nums == msg->knot_times.size());
            assert(num_knot_nums == msg->contact_surface_ids.size());

            for (size_t i = 0; i < msg->contact_surface_ids.size(); ++i)
            {
                std::vector<environment::SurfaceID> contact_surfaces_at_knot_i = msg->contact_surface_ids[i];

                assert(getEndEffectors().size() == contact_surfaces_at_knot_i.size());

                uint mask_bin = 0;
                for (size_t j = 0; j < contact_surfaces_at_knot_i.size(); ++j)
                {
                    mask_bin << 1;
                    if (contact_surfaces_at_knot_i[j] != galileo::legged::environment::NO_SURFACE)
                    {
                        mask_bin |= 1;
                    }
                }
            }

            setContactSequence(msg->knot_nums, msg->knot_times, mask_vec, msg->contact_surface_ids);
        }

        void GalileoLeggedRos::SurfaceCallback(const galileo_ros::msg::EnvironmentSurface::ConstPtr &msg)
        {
            // UNIMPLEMENTED
            addSurface(environment::createInfiniteGround());
        }

        void GalileoLeggedRos::InitializationCallback(const galileo_ros::msg::GalileoCommand::ConstPtr &msg)
        {
            auto X0 = casadi::DM(msg->X_initial);
            auto Xf = casadi::DM(msg->X_goal);

            assert(X0.size() == states()->nx);
            assert(Xf.size() == states()->nx);

            Initialize(X0, Xf);
        }

        void GalileoLeggedRos::UpdateCallback(const galileo_ros::msg::GalileoCommand::ConstPtr &msg)
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