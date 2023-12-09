#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace acro
{
    namespace environment
    {
        struct SurfaceData
        {
            // Assume it is Z aligned, so, all we need is the height.
            double origin_z_offset;

            // The A and b values defining the surface 2d polytope in the ground frame.
            Eigen::MatrixXd A;
            // A must be an "a by 2" vector.
            Eigen::VectorXd b;

            // The chebychev center of A and b
            Eigen::Vector2d polytope_local_chebyshev_center;
        };


        SurfaceData CreateInfiniteGround();

        template <class T>
        void PointViolation(const SurfaceData &region, const Eigen::Matrix<T, 2, 1>&point, Eigen::Matrix<T,Eigen::Dynamic,1> &ineq_violation);

        template <class T>
        void PointViolation(const SurfaceData &region, const Eigen::Matrix<T, 3, 1> &point, Eigen::Matrix<T,Eigen::Dynamic,1> &ineq_violation, Eigen::Matrix<T,Eigen::Dynamic,1> &eq_violation);

        bool isInRegion(const SurfaceData &region, const Eigen::Vector2d &point);

        bool isOnRegion(const SurfaceData &region, const Eigen::Vector3d &point);

        // The chebyshev center of A and b in 3d. WILL NOT WORK IF THE LOCAL CENTER HAS NOT BEEN COMPUTED.
        Eigen::Vector3d getChebyshevCenter(const SurfaceData &surface_data);

        Eigen::VectorXd CalculateChebyshevCenter(const Eigen::MatrixXd &A, const Eigen::VectorXd &b);

        using SurfaceID = int;

        const SurfaceID NO_SURFACE = -1;

        class EnvironmentSurfaces : public std::vector<SurfaceData>
        {
        public:
            EnvironmentSurfaces() : std::vector<SurfaceData>() {}

            std::vector<SurfaceID> getSurfacesUnder(const Eigen::Vector2d &ee_pos) const;

            std::vector<SurfaceData> getSurfacesFromIDs(const std::vector<SurfaceID> IDs) const;

            // Generate the straight shot trajectory of each limb from the starting to the target
            // and sample to find surfaces underneath

            // Get k-closest regions to current; convex program.
            std::vector<SurfaceID>
            getKClosestRegions(Eigen::Vector3d ee_pos, int k);

            uint num_surface() const { return (*this).size(); }
        };

    }
}