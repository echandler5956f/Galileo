#ifndef __galileo_multibody_force_base_hpp__
#define __galileo_multibody_force_base_hpp__

#include <pinocchio/multibody/data.hpp>
#include <pinocchio/spatial/force.hpp>

#include "galileo/multibody/fwd.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

#define GALILEO_FORCE_DATA_TYPEDEF(Force)                                  \
    using MatrixNcNv_t = typename traits<Force>::MatrixNcNv_t;             \
    using MatrixNcNdx_t = typename traits<Force>::MatrixNcNdx_t;           \
    using MatrixNcNu_t = typename traits<Force>::MatrixNcNu_t;

namespace galileo
{

    template <typename Derived, typename PhaseSpec>
    struct ForceDataBase : internal::CRTP<ForceDataBase<Derived, PhaseSpec>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        GALILEO_ROBOT_SPEC_MASTER_TYPEDEF(PS::RS);

        // The meta type of the derived class
        using Meta_t = typename traits<Derived>::Meta_t;
        // The data type of the meta type of the derived class
        using Data_t = typename traits<Meta_t>::Data_t;
        // The model type of the meta type of the derived class
        using Model_t = typename traits<Meta_t>::Model_t;

        GALILEO_FORCE_DATA_TYPEDEF(Meta_t);

        // Accessors required by ForceDataBase
        FORWARD_ACCESSOR(RobotData_t *, robot);
        FORWARD_ACCESSOR(FrameIndex_t, frame);
        FORWARD_ACCESSOR(ReferenceFrame_t, type);
        FORWARD_ACCESSOR(SE3_t, jMf);
        FORWARD_ACCESSOR(MatrixNcNv_t, Jc);
        FORWARD_ACCESSOR(Force_t, f);
        FORWARD_ACCESSOR(Force_t, fext);
        FORWARD_ACCESSOR(MatrixNcNdx_t, df_dx);
        FORWARD_ACCESSOR(MatrixNcNu_t, df_du);

        /**We have to override the CRTP derived() method to return the 
          derived object because ForceDataBase is multi-level CRTP**/

        /** Return reference to this as derived object */
        inline Derived &derived() & noexcept
        {
            return *static_cast<Derived *>(this);
        }
        /** Return reference to this as derived object */
        inline const Derived &derived() const & noexcept
        {
            return *static_cast<Derived const *>(this);
        }
        /** Return reference to this as derived object, when this is rvalue */
        inline Derived &&derived() && noexcept
        {
            return std::move(*static_cast<Derived *>(this));
        }

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

} // namespace galileo

#endif // __galileo_multibody_force_base_hpp__
