#include <casadi/casadi.hpp>
#include <string>
#include <iostream>
#include <fstream>


/**
 * @brief Saves the last solution to a file text in the form "t,x1,x2...,xn" where each row corresponds to a unique time step
 * 
 * @param last_solution The last solution from the solver
 * @param folder_location The location of the folder to save the file to
*/
void save_last_solution_to_file(casadi::MXVector& last_solution, std::string folder_location){
    std::ofstream file;
    file.open(folder_location + "/last_solution.txt");

    std::vector<double> times = traj.get_global_times();

    auto solx = last_solution[0];

    syd::vector<std::vector<double>> x_all;
    for(int i = 0; i < solx.rows(); i++){
        std::vector<double> xi_all = casadi::MX::evalf(solx(i,casadi::Slice(0, solx.columns()))).get_elements();
        x_all.push_back(xi_all);
    }

    for(int i = 0; i < x_all.size(); i++){
        file << times[i] << ", ";
        for(int j = 0; j < x_all[i].size()-1; j++){
            file << x_all[i][j] << ", ";
        }
        file << x_all[i][x_all[i].size()-1];
        file << std::endl;
    }   

    file.close();
};