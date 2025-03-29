#ifndef __galileo_core_activations_activation_model_base_hpp__
#define __galileo_core_activations_activation_model_base_hpp__

#include "galileo/core/activations/activation-base.hpp"

#define GALILEO_ROBOT_BASIC_TYPEDEF(Robot)                               \
    using Scalar = typename traits<Robot>::Scalar;                       \
    using VarScalar = typename traits<Robot>::VarScalar;                 \
    static constexpr int Options = traits<Robot>::Options;               \
    using RobotModelDerived = typename traits<Robot>::RobotModelDerived; \
    using RobotDataDerived = typename traits<Robot>::RobotDataDerived;

#define GALILEO_ROBOT_CONSTANTS(Robot)             \
    static constexpr int NX = traits<Robot>::NX;   \
    static constexpr int NU = traits<Robot>::NU;   \
    static constexpr int NDX = traits<Robot>::NDX; \
    static constexpr int NQ = traits<Robot>::NQ;   \
    static constexpr int NV = traits<Robot>::NV;

#define GALILEO_ROBOT_MODEL_TYPEDEF(Robot)

#define GALILEO_ROBOT_DATA_TYPEDEF(Robot)

namespace galileo
{
    namespace multibody
    {

        template <typename Derived>
        class RobotModelBase : internal::CRTP<Derived>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using RobotDerived = typename traits<Derived>::RobotDerived;
            GALILEO_ROBOT_BASIC_TYPEDEF(RobotDerived);
            GALILEO_ROBOT_CONSTANTS(RobotDerived);
            GALILEO_ROBOT_MODEL_TYPEDEF(RobotDerived);

        protected:
            inline RobotModelBase()
            {
            }

            inline RobotModelBase(const RobotModelBase &clone)
            {
                *this = clone;
            }

            inline RobotModelBase &operator=(const RobotModelBase &clone)
            {
                return *this;
            }

        }; // class RobotModelBase

    } // namespace multibody

} // namespace galileo

#endif // __galileo_multibody_robot_robot_model_base_hpp__
