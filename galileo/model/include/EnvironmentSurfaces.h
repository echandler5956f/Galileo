#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace acro
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
             * @brief Assume it is Z aligned, so, all we need is the height.
             *
             */
            double origin_z_offset;

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
             * @brief The chebychev center of A and b
             *
             */
            Eigen::Vector2d polytope_local_chebyshev_center;
        };

        /**
         * @brief Create a Infinite Ground object
         *
         * @return SurfaceData
         */
        SurfaceData CreateInfiniteGround();

        /**
         * @brief
         *
         * @tparam T
         * @param region
         * @param point
         * @param ineq_violation
         */
        template <class T>
        void PointViolation(const SurfaceData &region, const Eigen::Matrix<T, 2, 1> &point, Eigen::Matrix<T, Eigen::Dynamic, 1> &ineq_violation);

        /**
         * @brief
         *
         * @tparam T
         * @param region
         * @param point
         * @param ineq_violation
         * @param eq_violation
         */
        template <class T>
        void PointViolation(const SurfaceData &region, const Eigen::Matrix<T, 3, 1> &point, Eigen::Matrix<T, Eigen::Dynamic, 1> &ineq_violation, Eigen::Matrix<T, Eigen::Dynamic, 1> &eq_violation);

        /**
         * @brief
         *
         * @param region
         * @param point
         * @return true
         * @return false
         */
        bool isInRegion(const SurfaceData &region, const Eigen::Vector2d &point);

        /**
         * @brief
         *
         * @param region
         * @param point
         * @return true
         * @return false
         */
        bool isOnRegion(const SurfaceData &region, const Eigen::Vector3d &point);

        /**
         * @brief The chebyshev center of A and b in 3d. WILL NOT WORK IF THE LOCAL CENTER HAS NOT BEEN COMPUTED.
         *
         * @param surface_data
         * @return Eigen::Vector3d
         */
        Eigen::Vector3d getChebyshevCenter(const SurfaceData &surface_data);

        /**
         * @brief
         *
         * @param A
         * @param b
         * @return Eigen::VectorXd
         */
        Eigen::VectorXd CalculateChebyshevCenter(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);

        /**
         * @brief
         *
         */
        using SurfaceID = int;

        /**
         * @brief
         *
         */
        const SurfaceID NO_SURFACE = -1;

        /**
         * @brief
         *
         */
        class EnvironmentSurfaces : public std::vector<SurfaceData>
        {
        public:
            /**
             * @brief Construct a new Environment Surfaces object
             *
             */
            EnvironmentSurfaces() : std::vector<SurfaceData>() {}

            /**
             * @brief Generate the straight shot trajectory of each limb from the starting to the target and sample to find surfaces underneath.
             *
             * @param ee_pos
             * @return std::vector<SurfaceID>
             */
            std::vector<SurfaceID> getSurfacesUnder(const Eigen::Vector2d &ee_pos) const;

            /**
             * @brief Get the Surfaces From I Ds object
             *
             * @param IDs
             * @return std::vector<SurfaceData>
             */
            std::vector<SurfaceData> getSurfacesFromIDs(const std::vector<SurfaceID> IDs) const;

            /**
             * @brief Get k-closest regions to current; convex program.
             *
             * @param ee_pos
             * @param k
             * @return std::vector<SurfaceID>
             */
            std::vector<SurfaceID> getKClosestRegions(Eigen::Vector3d ee_pos, int k);

            /**
             * @brief
             *
             * @return uint
             */
            uint num_surface() const { return (*this).size(); }
        };

    }
}