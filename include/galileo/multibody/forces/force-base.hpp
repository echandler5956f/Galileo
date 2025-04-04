#ifndef __galileo_multibody_forces_force_base_hpp__
#define __galileo_multibody_forces_force_base_hpp__

#include <pinocchio/multibody/data.hpp>
#include <pinocchio/spatial/force.hpp>

#include "galileo/multibody/fwd.hpp"

namespace galileo
{
    namespace multibody
    {

        template <typename Derived, typename PhaseSpec>
        struct ForceDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            // RobotData_t *pinocchio;
            // pinocchio::FrameIndex frame;
            // pinocchio::ReferenceFrame type;
            // SE3 jMf;
            // MatrixNcNv_t Jc;
            // Force f;
            // Force fext;
            // MatrixNcNdx_t df_dx;
            // MatrixNcNu_t df_du;

        protected:
            inline ForceDataBase()
            {
            }

            inline ForceDataBase(const ForceDataBase &clone)
            {
                *this = clone;
            }

            inline ForceDataBase &operator=(const ForceDataBase &clone)
            {
                return *this;
            }

        }; // struct ForceDataBase

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_forces_force_base_hpp__
