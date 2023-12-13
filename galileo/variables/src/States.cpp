#include "galileo/variables/States.h"

namespace galileo
{
    namespace variables
    {
        States::States(const int nq_, const int nv_)
        {
            this->nq = nq_;
            this->nv = nv_;
            this->nx = this->nh + this->ndh + this->nq + this->nv;
            this->ndx = this->nh + this->ndh + 2 * this->nv;
            this->nvju = this->nv - this->nvb;
            this->nu = this->nF + this->nvju;
        }
    }
}