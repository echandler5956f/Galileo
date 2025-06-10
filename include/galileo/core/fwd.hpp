#ifndef __galileo_core_fwd_hpp__
#define __galileo_core_fwd_hpp__

#include "galileo/fwd.hpp"

namespace galileo
{

    enum AssignmentOp
    {
        setto,
        addto,
        rmfrom
    }; // enum AssignmentOp

    inline bool is_a_AssignmentOp(AssignmentOp op)
    {
        return (op == setto || op == addto || op == rmfrom);
    }

} // namespace galileo

#endif // __galileo_core_fwd_hpp__