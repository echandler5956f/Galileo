#ifndef __galileo_predictive_segments_segment_erk_hpp__
#define __galileo_predictive_segments_segment_erk_hpp__

#include "galileo/predictive/segments/segment-base.hpp"

namespace galileo
{
    namespace predictive
    {

        template <typename VarScalar,
                  typename NumScalar,
                  int Options,
                  int Ndeg,
                  ERKType _ERKType,
                  template <typename V, typename N, int O, int Ndeg> class ControlParamModelTpl,
                  template <typename V, typename N, int O> class NodeModelTpl>
        struct SegmentERKTpl;

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _Ndeg,
                  ERKType _ERKType,
                  template <typename V, typename N, int O, int _Ndeg> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        struct traits<SegmentERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>>
        {
            using ControlParamDerived = traits<_ControlParamModelTpl<_VarScalar, _NumScalar, _Options>>::ControlParamDerived;
            using NodeDerived = traits<_NodeModelTpl<_VarScalar, _NumScalar, _Options>>::NodeDerived;

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NX = traits<NodeDerived>::NX;
            static constexpr int NU = traits<NodeDerived>::NU;
            static constexpr int NDX = traits<NodeDerived>::NDX;
            static constexpr int NStages = _ERKType;

            using ControlParamModelDerived = traits<ControlParamDerived>::ControlParamModelDerived;
            using ControlParamDataDerived = traits<ControlParamDerived>::ControlParamDataDerived;

            using StageCoefficients_t = Eigen::Matrix<NumScalar, NStages, NStages, Options>;
            using Quadrature_t = Eigen::Matrix<NumScalar, NStages, 1, Options>;
            using Timings_t = Eigen::Matrix<NumScalar, NStages, 1, Options>;

            using SegmentDataDerived = SegmentDataERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>;
            using SegmentModelDerived = SegmentModelERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _Ndeg,
                  ERKType _ERKType,
                  template <typename V, typename N, int O, int _Ndeg> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        struct traits<SegmentDataERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>>
        {
            using SegmentDerived = SegmentERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>;
            using VarScalar = traits<SegmentDerived>::VarScalar;
            using NumScalar = traits<SegmentDerived>::NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _Ndeg,
                  ERKType _ERKType,
                  template <typename V, typename N, int O, int _Ndeg> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        struct traits<SegmentModelERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>>
        {
            using SegmentDerived = SegmentERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>;
            using VarScalar = traits<SegmentDerived>::VarScalar;
            using NumScalar = traits<SegmentDerived>::NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _Ndeg,
                  ERKType _ERKType,
                  template <typename V, typename N, int O, int _Ndeg> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        struct SegmentDataERKTpl : public SegmentDataERKBase<SegmentDataERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using SegmentDerived = SegmentERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>;
            using NodeDerived = typename traits<SegmentDerived>::NodeDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_ERK_BASIC_TYPEDEF(SegmentDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_ERK_CONSTANTS(SegmentDerived);

            GALILEO_NODE_DATA_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_DATA_ERK_TYPEDEF(SegmentDerived);

        }; // struct SegmentDataERKTpl

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  int _Ndeg,
                  ERKType _ERKType,
                  template <typename V, typename N, int O, int _Ndeg> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        class SegmentModelERKTpl : public SegmentModelERKBase<SegmentModelERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using SegmentDerived = SegmentERKTpl<_VarScalar, _NumScalar, _Options, _Ndeg, _ERKType, _ControlParamModelTpl, _NodeModelTpl>;
            using NodeDerived = typename traits<SegmentDerived>::NodeDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_ERK_BASIC_TYPEDEF(SegmentDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_ERK_CONSTANTS(SegmentDerived);

            GALILEO_NODE_MODEL_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_MODEL_ERK_TYPEDEF(SegmentDerived);

            template <typename StateVectorType, typename ControlMatrixType>
            void calc(SegmentDataDerived &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
            }

            template <typename StateVectorType, typename ControlMatrixType>
            void calcDiff(SegmentDataDerived &data,
                          const Eigen::MatrixBase<StateVectorType> &x,
                          const Eigen::MatrixBase<ControlMatrixType> &us) const
            {
            }

            template <typename StateVectorType, typename ControlMatrixType>
            void quasiStatic(SegmentDataDerived &data,
                             const Eigen::MatrixBase<StateVectorType> &x,
                             Eigen::MatrixBase<ControlMatrixType> &us,
                             const std::size_t maxiter, const Scalar tol) const
            {
            }

        protected:
            State_t state_;
            ControlParamModel_t control_param_model_;
            NodeModel_t node_model_;

            StageCoefficients_t stage_coefficients_;
            Quadrature_t quadrature_;
            Timings_t timings_;

        }; // class SegmentModelERKTpl

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_hpp__
