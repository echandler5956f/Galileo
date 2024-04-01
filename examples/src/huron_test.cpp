#include "huron_test.h"

int main(int argc, char **argv)
{
    std::vector<std::string> end_effector_names;
    std::vector<int> knot_num;
    std::vector<double> knot_time;
    std::vector<uint> mask_vec;
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;

    std::vector<double> q0_vec;
    std::vector<double> qf_vec;

    galileo::legged::helper::ReadProblemFromParameterFile(
        problem_parameter_location,
        end_effector_names,
        knot_num,
        knot_time,
        contact_surfaces,
        q0_vec,
        qf_vec);

    mask_vec = galileo::legged::helper::getMaskVectorFromContactSurfaces(contact_surfaces);

    galileo::legged::LeggedInterface solver_interface;

    solver_interface.LoadModel(huron_location, end_effector_names);
    solver_interface.LoadParameters(solver_parameter_location);

    casadi::DM X0;
    std::vector<double> X0_vec = galileo::legged::helper::getXfromq(solver_interface.states()->nx, solver_interface.states()->q_index, q0_vec);
    galileo::tools::vectorToCasadi(X0_vec, solver_interface.states()->nx, 1, X0);

    casadi::DM Xf;
    std::vector<double> Xf_vec = galileo::legged::helper::getXfromq(solver_interface.states()->nx, solver_interface.states()->q_index, qf_vec);
    galileo::tools::vectorToCasadi(Xf_vec, solver_interface.states()->nx, 1, Xf);

    solver_interface.addSurface(environment::createInfiniteGround());

    solver_interface.setContactSequence(
        knot_num,
        knot_time,
        mask_vec,
        contact_surfaces);

    solver_interface.Initialize(X0, Xf);
    solver_interface.Update(X0, Xf);

    Eigen::VectorXd new_times = Eigen::VectorXd::LinSpaced(250, 0., solver_interface.getRobotModel()->contact_sequence->getDT());
    Eigen::MatrixXd new_states = Eigen::MatrixXd::Zero(solver_interface.states()->nx, new_times.size());
    Eigen::MatrixXd new_inputs = Eigen::MatrixXd::Zero(solver_interface.states()->nu, new_times.size());
    // solver_interface.GetSolution(new_times, new_states, new_inputs);
    solver_interface.VisualizeSolutionAndConstraints(new_times, new_states, new_inputs);

    return 0;
}