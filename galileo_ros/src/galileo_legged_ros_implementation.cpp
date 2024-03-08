#include "galileo_ros/galileo_legged_ros_implementation.h"

// Constructor
GalileoLeggedROSImplementation::GalileoLeggedROSImplementation(ros::NodeHandle& node_handle, bool verbose)
    : node_handle_(node_handle), verbose_(verbose)
{
    // Initialize the subscribers
    InitSubscribers();

    // Initialize the publishers
    InitPublishers();

    // Initialize the services
    InitServices();
}


// Initialize the subscribers

void GalileoLeggedROSImplementation::InitSubscribers()
{
    // Subscribe to the robot model
    robot_model_subscriber_ = 
        node_handle_.subscribe("legged_robot_model", 1, &GalileoLeggedROSImplementation::LoadModelCallback, this);

    if(verbose_) ROS_INFO("Subscribed to legged_robot_model");

    // Subscribe to the parameter location strings
    parameter_location_subscriber_ = 
        node_handle_.subscribe("legged_parameter_location", 1, &GalileoLeggedROSImplementation::ParameterCallback, this);

    if(verbose_) ROS_INFO("Subscribed to legged_parameter_location");
    // // Subscribe to the robot state
    // robot_state_subscriber_ = 
    //     nh_.subscribe("legged_robot_state", 100, &GalileoLeggedROSImplementation::UpdateRobotState, this);

    // Subscribe to the robot command. this should end up being a service.
    // robot_command_subscriber_ = 
    //     nh_.subscribe("legged_robot_command", 100, &GalileoLeggedROSImplementation::UpdateRobotCommand, this);
}


// Initialize the publishers
void GalileoLeggedROSImplementation::InitPublishers()
{
    // Publish the solution 
    solution_publisher_ = 
        node_handle_.advertise<galileo_ros::RobotSolution>("legged_robot_solution", 100);
}

void GalileoLeggedROSImplementation::LoadModelCallback(const galileo_ros::ModelLocation::ConstPtr& model_location)
{
    // Load the robot model from the given model file and set the end effectors
    LoadModel(model_location->model_file_location, model_location->end_effector_names);

    robot_model_subscriber_.shutdown();
}

void GalileoLeggedROSImplementation::LoadModel(const std::string& model_file, const std::vector<std::string>& end_effector_names)
{
    std::string* end_effector_names_array = new std::string[end_effector_names.size()];

    for (size_t i = 0; i < end_effector_names.size(); ++i)
    {
        end_effector_names_array[i] = end_effector_names[i];
        if(verbose_) ROS_INFO("Setting end effector: %s", end_effector_names_array[i].c_str() );
    }
    if(verbose_) ROS_INFO ("Loading model from %s", model_file.c_str());
    
    // Load the robot model from the given model file and set the end effectors
    robot_ = std::make_shared<galileo::legged::LeggedBody>(model_file, end_effector_names.size(), end_effector_names_array);

    states_ = robot_->si;
}

void GalileoLeggedROSImplementation::ParameterCallback(const std_msgs::String::ConstPtr& msg)
{
    // RUN ONLY IF ROBOT MODEL IS LOADED
    if(!robot_)
    {
        ROS_ERROR("Robot model not loaded yet. Cannot update parameters.");
        return;
    }

    // Update the parameters from the given parameter file, and create the costs
    LoadParameters(msg->data);


    // get target position 
    // get cost weights
    // create cost function, L and Phi
    
    parameter_location_subscriber_.shutdown();
}


void GalileoLeggedROSImplementation::LoadParameters(const std::string& parameter_file)
{
    // opts["ipopt.linear_solver"] = "ma97";
    // opts["ipopt.ma97_order"] = "metis";
    opts_["ipopt.fixed_variable_treatment"] = "make_constraint";
    opts_["ipopt.max_iter"] = 250;


    if(verbose_) ROS_INFO("Parameters loaded from %s", parameter_file.c_str());
    // Do nothing, we will hard code parameters for now.
}



void GalileoLeggedROSImplementation::CreateProblemData()
{
    // First, get the costs L (running) and Phi (terminal). This creates the general problem data
    casadi::Function L, Phi;
    CreateCost(L, Phi);

    std::shared_ptr<GeneralProblemData> gp_data = std::make_shared<GeneralProblemData>(robot_->fint, robot_->fdif, L, Phi);

    // Init with a default gait and an infinite ground
    surfaces_ = std::make_shared<galileo::legged::environment::EnvironmentSurfaces>();
    surfaces_->push_back(galileo::legged::environment::createInfiniteGround());


    std::vector<int> knot_num = {20, 20, 20, 20, 20, 20, 20, 20, 20};
    std::vector<double> knot_time = {0.05, 0.3, 0.05, 0.3, 0.05, 0.3, 0.05, 0.3, 0.05};
    std::vector<uint> mask_vec = {0b1111, 0b1101, 0b1111, 0b1011, 0b1111, 0b1110, 0b1111, 0b0111, 0b1111}; // static walk

    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces = 
        {   {0, 0, 0, 0}, 
            {0, 0, galileo::legged::environment::NO_SURFACE, 0}, {0, 0, 0, 0}, 
            {0, galileo::legged::environment::NO_SURFACE, 0, 0}, {0, 0, 0, 0},
            {0, 0, 0, galileo::legged::environment::NO_SURFACE}, {0, 0, 0, 0}, 
            {galileo::legged::environment::NO_SURFACE, 0, 0, 0}, {0, 0, 0, 0},
            }; 


    int num_steps = 1;
    for (int i = 0; i < num_steps; ++i)
    {
        for (size_t j = 0; j < mask_vec.size(); ++j)
        {
            galileo::legged::contact::ContactMode mode;
            mode.combination_definition = robot_->getContactCombination(mask_vec[j]);
            mode.contact_surfaces = contact_surfaces[j];
            robot_->contact_sequence->addPhase(mode, knot_num[j], knot_time[j]);
        }
    }

    robot_->fillModeDynamics(false);

    // Create the problem data from the loaded model/parameters
    problem_data_ = std::make_shared<LeggedRobotProblemData>(gp_data,
                                                             surfaces_,
                                                             robot_->contact_sequence,
                                                             states_, std::make_shared<galileo::opt::ADModel>(robot_->cmodel),
                                                             std::make_shared<galileo::opt::ADData>(robot_->cdata), 
                                                             robot_->getEndEffectors(), 
                                                             robot_->cx, robot_->cu, robot_->cdt);
}

void GalileoLeggedROSImplementation::CreateCost( casadi::Function &L, casadi::Function &Phi ) const
{
    //HARDCODE INITIAL CONDITION 
    double q0[] = {
        0., 0., 0.339, 0., 0., 0., 1.,
        0., 0.67, -1.30, 0., 0.67, -1.3, 0., 0.67, -1.3, 0., 0.67, -1.3};

    Eigen::Map<galileo::opt::ConfigVector> q0_vec(q0, 19);

    casadi::DM X0 = casadi::DM::zeros(states_->nx, 1);
    int j = 0;
    for (int i = states_->nh + states_->ndh; i < states_->nh + states_->ndh + states_->nq; ++i)
    {
        X0(i) = q0_vec[j];
        ++j;
    }

    //HARDCODE COST WEIGHTS

    Eigen::VectorXd Q_diag(states_->nx); 
    Q_diag << 15., 15., 30., 5., 10., 10.,                          /*Centroidal momentum error weights*/
        0., 0., 0., 0., 0., 0.,                                     /*Rate of Centroidal momentum error weights*/
        500., 500., 500., 0.1, 0.1, 0.1,                            /*Floating base position and orientation (exponential coordinates) error weights*/
        20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., /*Joint position error weights*/
        0., 0., 0., 0., 0., 0.,                                     /*Floating base velocity error weights*/
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;             /*Joint velocity error weights*/
    Eigen::MatrixXd Q_mat = Q_diag.asDiagonal();

    Eigen::VectorXd R_diag(states_->nu);
    R_diag << 1e-3, 1e-3, 1e-3,                                     /*First contact wrench error weights*/
        1e-3, 1e-3, 1e-3,                                           /*Second contact wrench error weights*/
        1e-3, 1e-3, 1e-3,                                           /*Third contact wrench error weights*/
        1e-3, 1e-3, 1e-3,                                           /*Fourth contact wrench error weights*/
        10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10.; /*Joint acceleration error weights*/
    Eigen::MatrixXd R_mat = R_diag.asDiagonal();

    casadi::SX Q = casadi::SX::zeros(states_->ndx, states_->ndx);
    casadi::SX R = casadi::SX::zeros(states_->nu, states_->nu);

    pinocchio::casadi::copy(Q_mat, Q);
    pinocchio::casadi::copy(R_mat, R);

    //
    casadi::SX target_pos = casadi::SX::vertcat(casadi::SXVector{q0[0] + 0., q0[1] + 0., q0[2] + 0.});
    casadi::SX target_rot = casadi::SX::eye(3);

    pinocchio::SE3Tpl<galileo::opt::ADScalar, 0> oMf = robot_->cdata.oMf[robot_->model.getFrameId("base", pinocchio::BODY)];
    auto rot = oMf.rotation();
    auto pos = oMf.translation();
    casadi::SX crot = casadi::SX::zeros(3, 3);
    casadi::SX cpos = casadi::SX::zeros(3, 1);
    pinocchio::casadi::copy(rot, crot);
    pinocchio::casadi::copy(pos, cpos);

    casadi::SX rot_c = casadi::SX::mtimes(crot.T(), target_rot);

    casadi::SX target_error_casadi = casadi::SX::vertcat(casadi::SXVector{cpos - target_pos, casadi::SX::inv_skew(rot_c - rot_c.T()) / 2});

    casadi::SX X_ref = casadi::SX(X0);
    X_ref(casadi::Slice(states_->nh + states_->ndh, states_->nh + states_->ndh + 3)) = target_pos;
    casadi::SX X_error = casadi::SX::vertcat(casadi::SXVector{robot_->cx(casadi::Slice(0, states_->nh + states_->ndh)) - X_ref(casadi::Slice(0, states_->nh + states_->ndh)),
                                                  target_error_casadi,
                                                  robot_->cx(casadi::Slice(states_->nh + states_->ndh + states_->nqb, states_->nx)) - X_ref(casadi::Slice(states_->nh + states_->ndh + states_->nqb, states_->nx))});

    pinocchio::computeTotalMass(robot_->model, robot_->data);

    casadi::SX U_ref = casadi::SX::zeros(states_->nu, 1);
    casadi::SX u_error = robot_->cu - U_ref;

    L = casadi::Function("L",
                       {robot_->cx, robot_->cu},
                       {1. * 0.5 * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error)) +
                        1. * 0.5 * casadi::SX::dot(u_error, casadi::SX::mtimes(R, u_error))});

    Phi = casadi::Function("Phi",
                         {robot_->cx},
                         {1. * casadi::SX::dot(X_error, casadi::SX::mtimes(Q, X_error))});

}

std::vector< GalileoLeggedROSImplementation::LeggedConstraintBuilderType > 
GalileoLeggedROSImplementation::getLeggedConstraintBuilders() const
{
    
    LeggedConstraintBuilderType friction_cone_constraint_builder =
        std::make_shared<galileo::legged::constraints::FrictionConeConstraintBuilder<LeggedRobotProblemData>>();

    LeggedConstraintBuilderType velocity_constraint_builder =
        std::make_shared<galileo::legged::constraints::VelocityConstraintBuilder<LeggedRobotProblemData>>();

    LeggedConstraintBuilderType contact_constraint_builder =
        std::make_shared<galileo::legged::constraints::ContactConstraintBuilder<LeggedRobotProblemData>>();
        
    return {velocity_constraint_builder, friction_cone_constraint_builder};

}


void GalileoLeggedROSImplementation::CreateTrajOptSolver(){
    if(!problem_data_){
        ROS_ERROR("Problem data not created yet. Cannot create trajectory optimizer.");
        return;
    }
    if(!robot_){
        ROS_ERROR("robot_ not created yet. Cannot create trajectory optimizer.");
        return;
    }
    
    auto constraint_builders = getLeggedConstraintBuilders();
    decision_builder_ = std::make_shared<galileo::legged::constraints::LeggedDecisionDataBuilder<LeggedRobotProblemData>>();

    trajectory_opt_ = std::make_shared<LeggedTrajOpt>(problem_data_, robot_->contact_sequence, constraint_builders, decision_builder_, opts_);
}

/**
 * @brief Updates the solution based on the initial state.
 * 
 * @param initial_state The initial state of the robot.
 */
void GalileoLeggedROSImplementation::UpdateSolution( T_ROBOT_STATE initial_state ){
    if(!trajectory_opt_){
        ROS_ERROR("Trajectory optimizer not created yet. Cannot update solution.");
        return;
    }

    // Update the initial state
    trajectory_opt_->initFiniteElements(1, initial_state);

    //@todo set the initial guess and solve
}

int main(int argc, char** argv)
{
    // Initialize the ROS node
    ros::init(argc, argv, "gallileo_legged_ros_implementation_node");

    // Create a ROS node handle
    ros::NodeHandle nh;

    GalileoLeggedROSImplementation galileo_legged_imp(nh, true);

    // Spin the ROS node
    ros::spin();

    return 0;
}