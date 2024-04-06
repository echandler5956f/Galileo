#include "galileo/legged-model/LeggedModelHelpers.h"

namespace
{
    std::vector<std::vector<galileo::legged::environment::SurfaceID>> readContactSurfaces(std::vector<std::string> contact_surfaces_str)
    {
        std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces;
        std::vector<galileo::legged::environment::SurfaceID> current_contact_surface_combination;

        for (int i = 0; i < contact_surfaces_str.size(); i++)
        {
            bool some_vector_read = false;
            if (contact_surfaces_str[i].find_first_of('(') != std::string::npos)
            {
                int parenthesis_index = contact_surfaces_str[i].find_first_of('(');
                contact_surfaces_str[i] = contact_surfaces_str[i].substr(parenthesis_index + 1);
            }

            if (contact_surfaces_str[i].find_last_of(')') != std::string::npos)
            {
                some_vector_read = true;
                int parenthesis_index = contact_surfaces_str[i].find_last_of(')');
                contact_surfaces_str[i] = contact_surfaces_str[i].substr(0, parenthesis_index);
            }

            current_contact_surface_combination.push_back(std::stoi(contact_surfaces_str[i]));

            if (some_vector_read)
            {
                contact_surfaces.push_back(current_contact_surface_combination);
                current_contact_surface_combination.clear();
            }
        }

        return contact_surfaces;
    }
}
namespace galileo
{
    namespace legged
    {
        namespace helper
        {

            std::vector<double> getXfromq(int nx, int q_index, std::vector<double> q)
            {
                std::cout << "nx: " << nx << std::endl;
                std::cout << "q_index: " << q_index << std::endl;
                std::cout << "q.size(): " << q.size() << std::endl;
                assert(q_index + q.size() <= nx);

                std::vector<double> X(nx, 0); // Initialize X with size nx, fill with 0s
                std::copy(q.begin(), q.end(), X.begin() + q_index);

                return X;
            }

            void ReadProblemFromParameterFile(std::string problem_parameter_file_name,
                                              std::vector<std::string> &end_effector_names,
                                              std::vector<int> &knot_num,
                                              std::vector<double> &knot_time,
                                              std::vector<std::vector<galileo::legged::environment::SurfaceID>> &contact_surfaces,
                                              std::vector<double> &q0_vec,
                                              std::vector<double> &qf_vec)
            {

                std::map<std::string, std::string> data_map = galileo::tools::readFromFile(problem_parameter_file_name);

                end_effector_names = galileo::tools::readAsVector(data_map["end_effector_names"]);
                assert(end_effector_names.size() > 0);

                std::vector<std::string> knot_num_str = galileo::tools::readAsVector(data_map["knot_num"]);
                std::transform(knot_num_str.begin(), knot_num_str.end(), std::back_inserter(knot_num), [](const std::string &str)
                               { return std::stoi(str); });
                assert(knot_num.size() > 0);

                std::vector<std::string> knot_time_str = galileo::tools::readAsVector(data_map["knot_time"]);
                std::transform(knot_time_str.begin(), knot_time_str.end(), std::back_inserter(knot_time), [](const std::string &str)
                               { return std::stod(str); });
                assert(knot_time.size() == knot_num.size());

                std::vector<std::string> contact_surfaces_str = galileo::tools::readAsVector(data_map["contact_combinations"]);

                contact_surfaces = readContactSurfaces(contact_surfaces_str);
                assert(contact_surfaces.size() == knot_num.size());

                std::vector<std::string> q0_str = galileo::tools::readAsVector(data_map["q0"]);
                std::vector<std::string> qf_str = galileo::tools::readAsVector(data_map["qf"]);
                assert(q0_str.size() > 0);
                assert(q0_str.size() == qf_str.size());

                std::transform(q0_str.begin(), q0_str.end(), std::back_inserter(q0_vec), [](const std::string &str)
                               { return std::stod(str); });
                std::transform(qf_str.begin(), qf_str.end(), std::back_inserter(qf_vec), [](const std::string &str)
                               { return std::stod(str); });
            }

            std::vector<uint> getMaskVectorFromContactSurfaces(std::vector<std::vector<galileo::legged::environment::SurfaceID>> contact_surfaces)
            {
                std::vector<uint> mask_vec;
                int end_effector_size = contact_surfaces[0].size();
                for (auto &phase : contact_surfaces)
                {
                    assert(phase.size() == end_effector_size);

                    uint mask_bin = 0;
                    for (auto &csurface_id : phase)
                    {
                        mask_bin /= 2;
                        if (csurface_id != galileo::legged::environment::NO_SURFACE)
                        {
                            mask_bin += pow(2, end_effector_size - 1);
                        }
                    }

                    mask_vec.push_back(mask_bin);
                }
                return mask_vec;
            }
        }
    }
}
