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
                  template <typename V, typename N, int O> class ControlParamModelTpl,
                  template <typename V, typename N, int O> class NodeModelTpl>
        struct SegmentERKTpl;

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        struct traits<SegmentERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>>
        {
            using ControlParamDerived = traits<_ControlParamModelTpl<_VarScalar, _NumScalar, _Options>>::ControlParamDerived;
            using NodeDerived = traits<_NodeModelTpl<_VarScalar, _NumScalar, _Options>>::NodeDerived;

            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            static constexpr int NR = traits<ControlParamDerived>::NR;

            using SegmentDataDerived = SegmentDataERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>;
            using SegmentModelDerived = SegmentModelERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        struct traits<SegmentDataERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>>
        {
            using SegmentDerived = SegmentERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>;
            using VarScalar = traits<SegmentDerived>::VarScalar;
            using NumScalar = traits<SegmentDerived>::NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        struct traits<SegmentModelERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>>
        {
            using SegmentDerived = SegmentERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>;
            using VarScalar = traits<SegmentDerived>::VarScalar;
            using NumScalar = traits<SegmentDerived>::NumScalar;
        };

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        class SegmentModelERKTpl : public SegmentModelERKBase<SegmentModelERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using SegmentDerived = SegmentERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>;
            using NodeDerived = typename traits<SegmentDerived>::NodeDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_ERK_BASIC_TYPEDEF(SegmentDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_ERK_CONSTANTS(SegmentDerived);

            GALILEO_NODE_MODEL_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_MODEL_ERK_TYPEDEF(SegmentDerived);

        }; // class SegmentModelERKTpl

        template <typename _VarScalar,
                  typename _NumScalar,
                  int _Options,
                  template <typename V, typename N, int O> class _ControlParamModelTpl,
                  template <typename V, typename N, int O> class _NodeModelTpl>
        struct SegmentDataERKTpl : public SegmentDataERKBase<SegmentDataERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using SegmentDerived = SegmentERKTpl<_VarScalar, _NumScalar, _Options, _ControlParamModelTpl, _NodeModelTpl>;
            using NodeDerived = typename traits<SegmentDerived>::NodeDerived;

            GALILEO_NODE_BASIC_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_ERK_BASIC_TYPEDEF(SegmentDerived);

            GALILEO_NODE_CONSTANTS(NodeDerived);
            GALILEO_SEGMENT_ERK_CONSTANTS(SegmentDerived);

            GALILEO_NODE_DATA_TYPEDEF(NodeDerived);
            GALILEO_SEGMENT_DATA_ERK_TYPEDEF(SegmentDerived);

        }; // struct SegmentDataERKTpl

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_segments_segment_erk_hpp__
