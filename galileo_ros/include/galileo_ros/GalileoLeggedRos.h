#pragma once

#include <galileo/legged-model/LeggedInterface.h>

#include <ros/ros.h>
#include <ros/package.h>
#include <std_msgs/String.h>
#include <std_msgs/Float64.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_srvs/Trigger.h>

#include "galileo_ros/RobotModel.h"
#include "galileo_ros/ParameterFileLocation.h"
#include "galileo_ros/ContactSequence.h"
#include "galileo_ros/EnvironmentSurface.h"
#include "galileo_ros/GalileoCommand.h"
#include "galileo_ros/SolutionRequest.h"

namespace galileo
{
    namespace legged
    {
        class GalileoLeggedRos : public LeggedInterface
        {
        public:
            GalileoLeggedRos(std::shared_ptr<ros::NodeHandle> nh, std::string solver_id) : LeggedInterface(), nh_(nh), solver_id_(solver_id)
            {
                InitPublishers();
                InitSubscribers();
            }

            ~GalileoLeggedRos() {}

        private:
            // Initialize ROS

            // Initialize ROS publishers
            void InitPublishers();
            // Initialize ROS subscribers
            void InitSubscribers();
            // Initialize ROS services
            void InitServices();

            // Callbacks
            // Callback for model location subscriber, loads the model and sets the end effectors
            void ModelLocationCallback(const galileo_ros::RobotModel::ConstPtr &msg);
            // Callback for parameter location subscriber, loads the parameters
            void ParameterLocationCallback(const galileo_ros::ParameterFileLocation::ConstPtr &msg);
            // Callback for contact sequence subscriber, sets the contact sequence
            void ContactSequenceCallback(const galileo_ros::ContactSequence::ConstPtr &msg);
            // Callback for environment surface subscriber, adds an environment surface to the library of surfaces
            void SurfaceCallback(const galileo_ros::EnvironmentSurface::ConstPtr &msg);

            // Callback for initialization service, checks if the solver can be initialized
            bool CanInitServiceCallback(std_srvs::Trigger::Request &req, std_srvs::Trigger::Response &res);

            void GeneralCommandCallback(const galileo_ros::GalileoCommand::ConstPtr &msg)
            {
                if (msg->command_type == "init")
                {
                    InitializationCallback(msg);
                }
                UpdateCallback(msg);
            }

            // Callback for initialization command subscriber, initializes the solver
            void InitializationCallback(const galileo_ros::GalileoCommand::ConstPtr &msg);
            // Callback for update command subscriber, updates the solver
            void UpdateCallback(const galileo_ros::GalileoCommand::ConstPtr &msg);

            bool GetSolutionCallback(galileo_ros::SolutionRequest::Request &req, galileo_ros::SolutionRequest::Response &res);

            std::shared_ptr<ros::NodeHandle> nh_;
            std::string solver_id_;

            ros::Subscriber model_location_subscriber_;
            ros::Subscriber parameter_location_subscriber_;
            ros::Subscriber contact_sequence_subscriber_;
            ros::Subscriber surface_subscriber_;

            ros::Subscriber command_subscriber_;

            ros::ServiceServer can_init_service_;
            ros::ServiceServer get_solution_service_;
        };

    }
}
