#ifndef __galileo_multibody_end_effectors_hpp__
#define __galileo_multibody_end_effectors_hpp__

#include "galileo/multibody/fwd.hpp"

namespace galileo
{

    struct EndEffector
    {
        // Name of the end-effector frame
        std::string frame_name;

        // Index of the frame in the multibody model
        std::size_t frame_idx;

    }; // struct EndEffector

} // namespace galileo

#endif // __galileo_multibody_end_effectors_hpp__
