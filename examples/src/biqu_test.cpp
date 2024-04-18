#include "biqu_test.h"


std::vector<Eigen::Vector3d> loadSurfacePoints(std::string file_path)
{
    std::vector<Eigen::Vector3d> points;
    std::ifstream file;

    file.open(file_path);
    if (!file.is_open())
    {
        return points;
    }

    std::string line;
    std::cout << "Reading file " << file_path << std::endl;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        double x, y, z;
        if (!(iss >> x >> y >> z))
        {
            continue;
        }
        std::cout << "Read point " << x << " " << y << " " << z << std::endl;
        Eigen::Vector3d point;
        point[0] = x;
        point[1] = y;
        point[2] = z;
        points.push_back( point );
    }

    file.close();

    return points;
}

int main(int argc, char **argv)
{
    std::vector<std::string> end_effector_names;
    std::vector<int> knot_num;
    std::vector<double> knot_time;
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;

    std::vector<double> q0_vec;
    std::vector<double> qf_vec;

    galileo::legged::helper::ReadProblemFromParameterFile(problem_parameter_location,
                                                          end_effector_names,
                                                          knot_num,
                                                          knot_time,
                                                          contact_surfaces,
                                                          q0_vec,
                                                          qf_vec);

    galileo::legged::LeggedInterface solver_interface;

    solver_interface.LoadModel(robot_location, end_effector_names);
    solver_interface.LoadParameters(solver_parameter_location);


    // Load the environment surfaces
    std::filesystem::recursive_directory_iterator it(surface_folder);
    std::vector<std::string> files;
    for (auto &entry : it)
    {
        if (entry.is_regular_file())
        {
            files.push_back(entry.path().filename());
        }
    }

    std::vector<Eigen::Vector3d> points;
    for (auto file : files)
    {
        std::string file_path = surface_folder +"/" + file;
        std::vector<Eigen::Vector3d> points = loadSurfacePoints(file_path);
        environment::SurfaceData surface;
        environment::CreateSurfaceFromPolygon(points, surface);
        solver_interface.addSurface( surface );
    }

    int nx = solver_interface.states()->nx;
    int q_idx = solver_interface.states()->q_index;

    casadi::DM X0;
    std::vector<double> X0_vec = galileo::legged::helper::getXfromq(solver_interface.states()->nx, q_idx, q0_vec);
    galileo::tools::vectorToCasadi(X0_vec, nx, 1, X0);

    casadi::DM Xf;
    std::vector<double> Xf_vec = galileo::legged::helper::getXfromq(solver_interface.states()->nx, q_idx, qf_vec);
    galileo::tools::vectorToCasadi(Xf_vec, nx, 1, Xf);

    std::vector<uint> mask_vec = galileo::legged::helper::getMaskVectorFromContactSurfaces(contact_surfaces);

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