#ifndef __galileo_core_actuations_floating-base_hpp__
#define __galileo_core_actuations_floating-base_hpp__

#include "galileo/core/actuations/actuation-model-base.hpp"

namespace galileo
{

    namespace core
    {

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  int NX,
                  int NU,
                  int NDX>
        struct ActuationFloatingBaseTpl;

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NX,
                  int _NU,
                  int _NDX>
        struct traits<ActuationFloatingBaseTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>>
        {

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NX = _NX;
            static constexpr int NU = _NU;
            static constexpr int NDX = _NDX;

            using VectorNX_t = Eigen::Matrix<VarScalar, NX, 1, Options>;
            using VectorNU_t = Eigen::Matrix<VarScalar, NU, 1, Options>;
            using VectorNDX_t = Eigen::Matrix<VarScalar, NDX, 1, Options>;
            using MatrixNDX_t = Eigen::Matrix<VarScalar, NDX, NDX, Options>;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _NX,
                  int _NU,
                  int _NDX>
        class ActuationFloatingBaseTpl : public ActuationModelBase<ActuationFloatingBaseTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using ActuationDerived = ActuationFloatingBaseTpl<_VarScalar, _NumScalar, _Options, _NX, _NU, _NDX>;
            GALILEO_ACTUATIONS_BASIC_TYPEDEF(ActuationDerived);
            GALILEO_ACTUATIONS_CONSTANTS(ActuationDerived);
            GALILEO_ACTUATIONS_TYPEDEF(ActuationDerived);


        }; // class ActuationFloatingBaseTpl

    } // namespace core

} // namespace galileo

#endif // __galileo_core_actuations_floating-base_hpp__