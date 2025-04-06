#ifndef __galileo_multibody_force_base_hpp__
#define __galileo_multibody_force_base_hpp__

#include "galileo/multibody/fwd.hpp"

#define GALILEO_FORCE_DATA_TYPEDEF(Force)                                  \
    using RobotDataPointer_t = typename traits<Force>::RobotDataPointer_t; \
    using Index_t = typename traits<Force>::Index_t;                       \
    using ReferenceFrame_t = typename traits<Force>::ReferenceFrame_t;     \
    using SE3_t = typename traits<Force>::SE3_t;                           \
    using MatrixNcNv_t = typename traits<Force>::MatrixNcNv_t;             \
    using Force_t = typename traits<Force>::Force_t;                       \
    using MatrixNcNdx_t = typename traits<Force>::MatrixNcNdx_t;           \
    using MatrixNcNu_t = typename traits<Force>::MatrixNcNu_t;

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
            using ForceDerived = typename traits<Derived>::ForceDerived;

            GALILEO_FORCE_DATA_TYPEDEF(ForceDerived);

            FORWARD_ACCESSOR(RobotDataPointer_t, robot_data_pointer);
            FORWARD_ACCESSOR(Index_t, frame);
            FORWARD_ACCESSOR(ReferenceFrame_t, type);
            FORWARD_ACCESSOR(SE3_t, jMf);
            FORWARD_ACCESSOR(MatrixNcNv_t, Jc);
            FORWARD_ACCESSOR(Force_t, f);
            FORWARD_ACCESSOR(Force_t, fext);
            FORWARD_ACCESSOR(MatrixNcNdx_t, df_dx);
            FORWARD_ACCESSOR(MatrixNcNu_t, df_du);

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

#endif // __galileo_multibody_force_base_hpp__
