#ifndef __galileo_multibody_spatial_force_base_hpp__
#define __galileo_multibody_spatial_force_base_hpp__

#include <pinocchio/multibody/data.hpp>
#include <pinocchio/spatial/force.hpp>

#include "galileo/domains/multibody/spatial/fwd.hpp"

#define GALILEO_FORCE_DATA_TYPEDEF(Force) \
    using MatrixNcNv_t = typename traits<Force>::MatrixNcNv_t; \
    using MatrixNcNdx_t = typename traits<Force>::MatrixNcNdx_t; \
    using MatrixNcNu_t = typename traits<Force>::MatrixNcNu_t;

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct ForceDataBase : public internal::CRTP<Derived>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        GALILEO_FORCE_DATA_TYPEDEF(Meta_t);

        // Accessors required by ForceDataBase
        FORWARD_ACCESSOR(typename PS::RobotData_t *, robot);
        FORWARD_ACCESSOR(typename PS::FrameIndex_t, frame);
        FORWARD_ACCESSOR(typename PS::ReferenceFrame_t, type);
        FORWARD_ACCESSOR(typename PS::SE3_t, jMf);
        FORWARD_ACCESSOR(typename PS::Force_t, f);
        FORWARD_ACCESSOR(typename PS::Force_t, fext);

        FORWARD_ACCESSOR(MatrixNcNv_t, Jc);
        FORWARD_ACCESSOR(MatrixNcNdx_t, df_dx);
        FORWARD_ACCESSOR(MatrixNcNu_t, df_du);

    protected:
        inline ForceDataBase() {}
        inline ForceDataBase(const ForceDataBase &clone) { *this = clone; }
        inline ForceDataBase &operator=(const ForceDataBase &clone) { return *this; }

    }; // struct ForceDataBase

} // namespace galileo

#endif // __galileo_multibody_spatial_force_base_hpp__
