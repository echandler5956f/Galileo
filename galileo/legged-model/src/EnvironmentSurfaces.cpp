#include "galileo/legged-model/EnvironmentSurfaces.h"

namespace galileo
{
    namespace legged
    {
        namespace environment
        {

            SurfaceData createInfiniteGround()
            {
                SurfaceData infinite_ground;
                infinite_ground.surface_transform = Eigen::Transform<double, 3, Eigen::Affine>::Identity();
                infinite_ground.A = {0, 0};

                infinite_ground.b = Eigen::VectorXd::Zero(1);
                infinite_ground.polytope_local_chebyshev_center = {0, 0};
                return infinite_ground;
            }

            template <class T>
            void pointViolation(const SurfaceData &region, const Eigen::Matrix<T, 2, 1> &point, Eigen::Matrix<T, Eigen::Dynamic, 1> &ineq_violation)
            {
                ineq_violation = region.A * point - region.b;
            }
            template <class T>
            void pointViolation(const SurfaceData &region, const Eigen::Matrix<T, 3, 1> &point, Eigen::Matrix<T, Eigen::Dynamic, 1> &ineq_violation, Eigen::Matrix<T, 1, 1> &eq_violation)
            {
                auto point_in_surface = region.WorldToSurface(point);
                ineq_violation = region.A * point_in_surface.head(2) - region.b;
                eq_violation = point_in_surface.tail(1);
            }

            bool isInRegion(const SurfaceData &region, const Eigen::Vector2d &point)
            {
                Eigen::VectorXd violation;
                pointViolation(region, point, violation);
                return (violation.array() <= Eigen::VectorXd::Zero(violation.size()).array()).all();
            }

            bool isOnRegion(const SurfaceData &region, const Eigen::Vector3d &point)
            {
                Eigen::VectorXd ineq_violation;
                Eigen::Matrix<double, 1, 1> eq_violation;
                pointViolation(region, point, ineq_violation, eq_violation);
                bool ineq_satisfied = (ineq_violation.array() <= Eigen::VectorXd::Zero(ineq_violation.size()).array()).all();
                bool eq_satisfied = (eq_violation[0] <= 1e-6) && (eq_violation[0] >= -1e-6);
                return ineq_satisfied && eq_satisfied;
            }

            // The chebyshev center of A and b in 3d. WILL NOT WORK IF THE LOCAL CENTER HAS NOT BEEN COMPUTED.
            Eigen::Vector3d getChebyshevCenter(const SurfaceData &surface_data)
            {
                // returns the chebychev center in the global frame.
                // (Adds a height component if the surface is gravity aligned and in the global frame)

                Eigen::Vector3d c_center;
                c_center.head(2) = surface_data.polytope_local_chebyshev_center;
                c_center.tail(1)[0] = 0;
                return surface_data.surface_transform * c_center;
            }

            // Eigen::VectorXd CalculateCenter(const Eigen::MatrixXd &A, const Eigen::VectorXd &b)
            // {
            //     // min (-r)
            //     // r (scalar)
            //     // d in Rn
            //     //  s.t
            //     //       (r * a_i / ||a_i||) + d  \leq  ||a_i|| * b
            //     // where a_i = A[i,:] .transpose()

            //     // or

            //     // min [[-1] [0 0 ... 0]] * x
            //     // x = [r;d]
            //     //  s.t [normalized(a_i) I_(nxn)] * r \leq  b~
            //     // where b~_i = b_i * ||a_i||

            //     // this is a linear program.
            // }

            std::vector<int> EnvironmentSurfaces::getSurfacesUnder(const Eigen::Vector2d &ee_pos) const
            {
                std::vector<int> surface_indeces;
                for (std::size_t i = 0; i < this->size(); i++)
                {
                    bool is_in_region = isInRegion((*this)[i], ee_pos);
                    if (is_in_region)
                    {
                        surface_indeces.push_back(i);
                    }
                }
                return surface_indeces;
            }

            SurfaceData EnvironmentSurfaces::getSurfaceFromID(const SurfaceID ID) const
            {
                return (*this)[ID];
            }

            std::vector<SurfaceData> EnvironmentSurfaces::getSurfacesFromIDs(const std::vector<int> IDs) const
            {
                std::vector<SurfaceData> surfaces;
                for (auto i : IDs)
                {
                    surfaces.push_back(getSurfaceFromID(i));
                }
                return surfaces;
            }

            // Generate the straight shot trajectory of each limb from the starting to the target
            // and sample to find surfaces underneath

            // Get k-closest regions to current; convex program.
            // std::vector<int>
            // LeggedBody::getKClosestRegions(Eigen::Vector3d ee_pos, int k) {}

        }
    }
}