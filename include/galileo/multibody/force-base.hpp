#ifndef __galileo_multibody_force_base_hpp__
#define __galileo_multibody_force_base_hpp__

#include "galileo/multibody/fwd.hpp"

#define GALILEO_FORCE_DATA_BASIC_TYPEDEF(ForceData)            \
    using VarScalar = typename traits<ForceData>::VarScalar;   \
    using NumScalar = typename traits<ForceData>::NumScalar;   \
    static constexpr int Options = traits<ForceData>::Options; \
    using ForceDataDerived = typename traits<ForceData>::ForceDataDerived;

#define GALILEO_FORCE_DATA_CONSTANTS(ForceData)        \
    static constexpr int NX = traits<ForceData>::NX;   \
    static constexpr int NU = traits<ForceData>::NU;   \
    static constexpr int NDX = traits<ForceData>::NDX; \
    static constexpr int NQ = traits<ForceData>::NQ;   \
    static constexpr int NV = traits<ForceData>::NV;   \
    static constexpr int NC = traits<ForceData>::NC;

#define GALILEO_FORCE_DATA_TYPEDEF(ForceData)                    \
    using RobotData_t = typename traits<ForceData>::RobotData_t; \
    using MatrixNCNV_t = typename traits<ForceData>::MatrixNCNV_t; \
    using MatrixNCNDX_t = typename traits<ForceData>::MatrixNCNDX_t; \
    using MatrixNCNU_t = typename traits<ForceData>::MatrixNCNU_t;

namespace galileo
{
    namespace multibody
    {

        template <typename Derived>
        struct ForceDataBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ForceDataDerived = typename traits<Derived>::ForceDataDerived;
            GALILEO_FORCE_DATA_BASIC_TYPEDEF(ForceDataDerived);
            GALILEO_FORCE_DATA_CONSTANTS(ForceDataDerived);
            GALILEO_FORCE_DATA_TYPEDEF(ForceDataDerived);

            // PinocchioData_t *pinocchio;
            // pinocchio::FrameIndex frame;
            // pinocchio::ReferenceFrame type;
            // SE3 jMf;
            // MatrixNCNV_t Jc;
            // Force f;
            // Force fext;
            // MatrixNCNDX_t df_dx;
            // MatrixNCNU_t df_du;

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
