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

// #define GALILEO_FORCE_DATA_BASE_DEFAULT_ACCESSOR(Force)      \
//     FORWARD_ACCESSOR(RobotDataPointer_t, robot_data_pointer) \
//     FORWARD_ACCESSOR(Index_t, frame)                         \
//     FORWARD_ACCESSOR(ReferenceFrame_t, type)                 \
//     FORWARD_ACCESSOR(SE3_t, jMf)                             \
//     FORWARD_ACCESSOR(MatrixNcNv_t, Jc)                       \
//     FORWARD_ACCESSOR(Force_t, f)                             \
//     FORWARD_ACCESSOR(Force_t, fext)                          \
//     FORWARD_ACCESSOR(MatrixNcNdx_t, df_dx)                   \
//     FORWARD_ACCESSOR(MatrixNcNu_t, df_du)

namespace galileo
{
    namespace multibody
    {

        template <typename Derived>
        struct ForceDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = typename traits<Derived>::PS;
            GALILEO_FORCE_DATA_TYPEDEF(Derived);

            FORWARD_ACCESSOR(RobotDataPointer_t, robot_data_pointer);
            FORWARD_ACCESSOR(Index_t, frame);
            FORWARD_ACCESSOR(ReferenceFrame_t, type);
            FORWARD_ACCESSOR(SE3_t, jMf);
            FORWARD_ACCESSOR(MatrixNcNv_t, Jc);
            FORWARD_ACCESSOR(Force_t, f);
            FORWARD_ACCESSOR(Force_t, fext);
            FORWARD_ACCESSOR(MatrixNcNdx_t, df_dx);
            FORWARD_ACCESSOR(MatrixNcNu_t, df_du);

            // const RobotDataPointer_t &robot_data() const
            // {
            //     return derived().robot_data_pointer_accessor();
            // }

            // RobotDataPointer_t &robot_data()
            // {
            //     return derived().robot_data_pointer_accessor();
            // }

            // const Index_t &frame() const
            // {
            //     return derived().frame_accessor();
            // }

            // Index_t &frame()
            // {
            //     return derived().frame_accessor();
            // }

            // const ReferenceFrame_t &type() const
            // {
            //     return derived().type_accessor();
            // }

            // ReferenceFrame_t &type()
            // {
            //     return derived().type_accessor();
            // }

            // const SE3_t &jMf() const
            // {
            //     return derived().jMf_accessor();
            // }

            // SE3_t &jMf()
            // {
            //     return derived().jMf_accessor();
            // }

            // const MatrixNcNv_t &Jc() const
            // {
            //     return derived().Jc_accessor();
            // }

            // MatrixNcNv_t &Jc()
            // {
            //     return derived().Jc_accessor();
            // }

            // const Force_t &f() const
            // {
            //     return derived().f_accessor();
            // }

            // Force_t &f()
            // {
            //     return derived().f_accessor();
            // }

            // const Force_t &fext() const
            // {
            //     return derived().fext_accessor();
            // }

            // Force_t &fext()
            // {
            //     return derived().fext_accessor();
            // }

            // const MatrixNcNdx_t &df_dx() const
            // {
            //     return derived().df_dx_accessor();
            // }

            // MatrixNcNdx_t &df_dx()
            // {
            //     return derived().df_dx_accessor();
            // }

            // const MatrixNcNu_t &df_du() const
            // {
            //     return derived().df_du_accessor();
            // }

            // MatrixNcNu_t &df_du()
            // {
            //     return derived().df_du_accessor();
            // }

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
