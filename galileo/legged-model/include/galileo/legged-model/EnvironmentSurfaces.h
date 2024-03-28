#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>

namespace galileo
{
    namespace legged
    {
        namespace environment
        {
            /**
             * @brief
             *
             */
            struct SurfaceData
            {
                /**
                 * @brief The transform of the surface in the global frame.
                 */
                Eigen::Transform<double, 3, Eigen::Affine> surface_transform;

                /**
                 * @brief The A and b values defining the surface 2d polytope in the ground frame.
                 *
                 */
                Eigen::MatrixXd A;

                /**
                 * @brief A must be an "a by 2" vector.
                 *
                 */
                Eigen::VectorXd b;

                /**
                 * @brief The chebychev center of A and b.
                 *
                 */
                Eigen::Vector2d polytope_local_chebyshev_center;

                /**
                 * @brief Get a "rotation" matrix R such that R * z_hat = unit_surface_normal_in_global_frame
                 */
                Eigen::Matrix<double, 3, 3> Rotation() const
                {
                    return surface_transform.rotation();
                }
                
                /**
                 * @brief Get the Inverse of the translation; used to transform points from the global frame to the surface frame.
                 * @param world_point The point to be transformed.
                 */
                Eigen::VectorXd WorldToSurface(Eigen::VectorXd world_point) const
                {
                    return Rotation().transpose() * (world_point - surface_transform.translation());
                }
            };

            /**
             * @brief Creates an infinite ground surface.
             *
             * This function creates an infinite ground surface using the SurfaceData class.
             *
             * @return The created infinite ground surface.
             */
            SurfaceData createInfiniteGround();

            /**
             * \brief Calculates the violation of inequality constraints at a given point within a surface region.
             *
             * This function calculates the violation of inequality constraints at a specified point within a surface region.
             * The violation is stored in the provided ineq_violation vector.
             *
             * \param region The surface region to evaluate.
             * \param point The point at which to calculate the violation.
             * \param ineq_violation The vector to store the violation of inequality constraints.
             *
             * \tparam T The data type of the matrix elements.
             */
            template <class T>
            void pointViolation(const SurfaceData &region, const Eigen::Matrix<T, 2, 1> &point, Eigen::Matrix<T, Eigen::Dynamic, 1> &ineq_violation);

            /**
             * @brief Calculates the violation of point constraints for a given region and point.
             *
             * This function calculates the violation of point constraints for a given region and point.
             * It updates the ineq_violation and eq_violation matrices with the calculated violations.
             *
             * @tparam T The type of the elements in the matrices.
             * @param region The surface data of the region.
             * @param point The point for which the violations are calculated.
             * @param ineq_violation The matrix to store the inequality violations.
             * @param eq_violation The matrix to store the equality violations.
             */
            template <class T>
            void pointViolation(const SurfaceData &region, const Eigen::Matrix<T, 3, 1> &point, Eigen::Matrix<T, Eigen::Dynamic, 1> &ineq_violation, Eigen::Matrix<T, Eigen::Dynamic, 1> &eq_violation);

            /**
             * @brief Checks if a given point is inside a specified region.
             *
             * @param region The surface region to check against.
             * @param point The point to check.
             * @return True if the point is inside the region, false otherwise.
             */
            bool isInRegion(const SurfaceData &region, const Eigen::Vector2d &point);

            /**
             * @brief Checks if a given point is on a specified region.
             *
             * @param region The surface data representing the region.
             * @param point The point to be checked.
             * @return True if the point is on the region, false otherwise.
             */
            bool isOnRegion(const SurfaceData &region, const Eigen::Vector3d &point);

            /**
             * @brief Calculates the Chebyshev center of a surface.
             *
             * This function calculates the Chebyshev center of a given surface using the provided surface data.
             *
             * @param surface_data The surface data used to calculate the Chebyshev center.
             * @return The Chebyshev center as a 3D vector.
             */
            Eigen::Vector3d getChebyshevCenter(const SurfaceData &surface_data);

            /**
             * @brief Calculates the Chebyshev center of a linear inequality system.
             *
             * This function calculates the Chebyshev center of a linear inequality system represented by the matrix A and vector b.
             * The Chebyshev center is the center of the largest ball that lies inside all the half-spaces defined by the inequalities.
             *
             * @param A The matrix representing the linear inequality system.
             * @param b The vector representing the right-hand side of the linear inequality system.
             * @return The Chebyshev center as a vector.
             */
            Eigen::VectorXd calculateChebyshevCenter(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);

            /**
             * @brief Represents the ID of a surface in the environment.
             */
            using SurfaceID = int;

            /**
             * @brief The ID of a surface that is not defined.
             *
             */
            const SurfaceID NO_SURFACE = -1;

            /**
             * @brief A vector of surface data.
             *
             */
            class EnvironmentSurfaces : public std::vector<SurfaceData>
            {
            public:
                /**
                 * @brief Construct a new Environment Surfaces object.
                 *
                 */
                EnvironmentSurfaces() : std::vector<SurfaceData>() {}

                /**
                 * @brief Generate the straight shot trajectory of each limb from the starting to the target and sample to find surfaces underneath.
                 *
                 * @param ee_pos The position of the end effector.
                 * @return std::vector<SurfaceID> The IDs of the surfaces underneath the end effector.
                 */
                std::vector<SurfaceID> getSurfacesUnder(const Eigen::Vector2d &ee_pos) const;

                /**
                 * @brief Get the surface data from the ID.
                 * 
                 * @param ID The ID of the surface.
                 * @return SurfaceData The surface data of the surface.
                 */
                SurfaceData getSurfaceFromID(const SurfaceID ID) const;

                /**
                 * @brief Get the surface data from the IDs.
                 *
                 * @param IDs The IDs of the surfaces.
                 * @return std::vector<SurfaceData> The surface data of the surfaces.
                 */
                std::vector<SurfaceData> getSurfacesFromIDs(const std::vector<SurfaceID> IDs) const;

                /**
                 * @brief Get k-closest regions to current; convex program.
                 *
                 * @param ee_pos The position of the end effector.
                 * @param k The number of closest regions to find.
                 * @return std::vector<SurfaceID> The IDs of the k-closest regions.
                 */
                std::vector<SurfaceID> getKClosestRegions(Eigen::Vector3d ee_pos, int k);

                /**
                 * @brief Get the number of surfaces.
                 *
                 * @return uint The number of surfaces.
                 */
                uint numSurface() const { return (*this).size(); }
            };

        }
    }
}