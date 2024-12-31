#pragma once

#include "galileo/core/node-base.hpp"

#define GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, TYPENAME)               \
    typedef TYPENAME traits<Segment>::Scalar Scalar;                           \
    typedef TYPENAME traits<Segment>::SegmentModelDerived SegmentModelDerived; \
    typedef TYPENAME traits<Segment>::SegmentDataDerived SegmentDataDerived;   \
    typedef TYPENAME traits<Segment>::NodeModel_t NodeModel_t;                 \
    typedef TYPENAME traits<Segment>::NodeData_t NodeData_t;                   \
    typedef TYPENAME traits<Segment>::StateModel_t StateModel_t;               \
    typedef TYPENAME traits<Segment>::ControlModel_t ControlModel_t;           \
    typedef TYPENAME traits<Segment>::VectorXs_t VectorXs_t;                   \
    typedef TYPENAME traits<Segment>::MatrixXs_t MatrixXs_t;

#define GALILEO_SEGMENT_TYPEDEF(Segment) GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, typename)
#define GALILEO_SEGMENT_TYPEDEF_TEMPLATE(Segment) \
    GALILEO_SEGMENT_MODEL_TYPEDEF_GENERIC(Segment, typename)

#define GALILEO_SEGMENT_CAST_TYPE_SPECIALIZATION(SegmentModelTpl) \
    template <typename Scalar, typename NewScalar>                \
    struct CastType<NewScalar, SegmentModelTpl<Scalar>>           \
    {                                                             \
        typedef SegmentModelTpl<NewScalar> type;                  \
    }

namespace galileo
{
    /**
     * @brief Abstract class for representing a segment of a trajectory
     *
     * A segment is defined by a list of nodes and a model that describes how
     * the nodes are connected. This segment model can either be shooting based,
     * where the nodes are connected in a sequential chain via explicit integration,
     * or collocation based, where the nodes are connected in a graph via implicit
     * integration.
     */
    template <typename Derived>
    class SegmentModelBase : NumericalBase<Derived>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::SegmentDerived SegmentDerived;
        GALILEO_SEGMENT_TYPEDEF_TEMPLATE(SegmentDerived);

        SegmentModelDerived &derived()
        {
            return *static_cast<Derived *>(this);
        }

        const SegmentModelDerived &derived() const
        {
            return *static_cast<const Derived *>(this);
        }

        void calc(const boost::shared_ptr<SegmentDataDerived> &data,
                  const Eigen::Ref<const MatrixXs_t> &x,
                  const Eigen::Ref<const MatrixXs_t> &u)
        {
            derived().calc(data, x, u);
        }

        void calc(const boost::shared_ptr<SegmentDataDerived> &data,
                  const Eigen::Ref<const MatrixXs_t> &x)
        {
            derived().calc(data, x);
        }

        void calcDiff(const boost::shared_ptr<SegmentDataDerived> &data,
                      const Eigen::Ref<const MatrixXs_t> &x,
                      const Eigen::Ref<const MatrixXs_t> &u)
        {
            derived().calcDiff(data, x, u);
        }

        void calcDiff(const boost::shared_ptr<SegmentDataDerived> &data,
                      const Eigen::Ref<const MatrixXs_t> &x)
        {
            derived().calcDiff(data, x);
        }

        void quasiStatic(const boost::shared_ptr<SegmentDataDerived> &data,
                         Eigen::Ref<MatrixXs_t> u, const Eigen::Ref<const MatrixXs_t> &x,
                         const std::size_t maxiter = 100, const Scalar tol = Scalar(1e-9))
        {
            derived().quasiStatic(data, u, x, maxiter, tol);
        }

        std::vector<VectorXs_t> quasiStatic_xs(const boost::shared_ptr<SegmentDataDerived> &data,
                                               const Eigen::Ref<const MatrixXs_t> &x,
                                               const std::size_t maxiter = 100,
                                               const Scalar tol = Scalar(1e-9))
        {
            return derived().quasiStatic_xs(data, x, maxiter, tol);
        }

        void circularAppend(boost::shared_ptr<NodeModel_t> model, boost::shared_ptr<NodeData_t> data)
        {
            derived().circularAppend(model, data);
        }

        void circularAppend(boost::shared_ptr<NodeModel_t> model)
        {
            derived().circularAppend(model);
        }

        void updateNode(const std::size_t &id, boost::shared_ptr<NodeModel_t> model, boost::shared_ptr<NodeData_t> data)
        {
            derived().updateNode(id, model, data);
        }

        void updateNode(const std::size_t &id, boost::shared_ptr<NodeModel_t> model)
        {
            derived().updateNode(id, model);
        }

        std::size_t getNumNodes() const
        {
            return derived().getNumNodes();
        }

        const std::vector<boost::shared_ptr<NodeModel_t>> &getNodes() const
        {
            return derived().getNodes();
        }

        const std::vector<boost::shared_ptr<NodeData_t>> &getNodesData() const
        {
            return derived().getNodesData();
        }

        const ControlModel_t &getControl() const
        {
            return derived().getControl();
        }

        const VectorXs_t &getNodeTimes() const
        {
            return derived().getNodeTimes();
        }

        Scalar getPeriod() const
        {
            return derived().getPeriod();
        }

        std::size_t getNx() const
        {
            return derived().getNx();
        }

        std::size_t getNdx() const
        {
            return derived().getNdx();
        }

        std::size_t getNu() const
        {
            return derived().getNu();
        }

        bool isUpdated() const
        {
            return derived().isUpdated();
        }

    protected:
        // Default constructor: protected
        // Prevent the construction of stand-alone SegmentModelBase.
        inline SegmentModelBase() : Ns_(0), h_(0), nx_(0), ndx_(0), nu_(0), is_updated_(false)
        {
        }

        // Copy constructor: protected
        // Copy of stand-alone SegmentModelBase are prevented, but can be used from inheriting
        // objects.
        inline SegmentModelBase(const SegmentModelBase &clone)
        {
            *this = clone;
        }

        // Copy operator: protected
        // Copy of stand-alone SegmentModelBase are prevented, but can be used from inheriting
        // objects.
        inline SegmentModelBase &operator=(const SegmentModelBase &clone)
        {
            Ns_ = clone.Ns_;
            nodes_ = clone.nodes_;
            nodes_data_ = clone.nodes_data_;
            control_ = clone.control_;
            node_times_ = clone.node_times_;
            h_ = clone.h_;
            nx_ = clone.nx_;
            ndx_ = clone.ndx_;
            nu_ = clone.nu_;
            is_updated_ = clone.is_updated_;
            return *this;
        }

        std::size_t Ns_;                                        //!< Number of nodes in this segment
        std::vector<boost::shared_ptr<NodeModel_t>> nodes_;     //!< List of nodes in this segment
        std::vector<boost::shared_ptr<NodeData_t>> nodes_data_; //!< List of nodes data in this segment
        ControlModel_t control_;                                //!< Control parameterization model for this segment
        VectorXs_t node_times_;                                 //!< Vector of node times, normalized to [0, 1]
        Scalar h_;                                              //!< Segment period

        std::size_t nx_;  //!< State dimension for this segment
        std::size_t ndx_; //!< State rate dimension for this segment
        std::size_t nu_;  //!< Control dimension for this segment
        bool is_updated_;
    };

    // A segment is either a shooting or collocation segment.
    // Shooting segments have block diagonal C matrices, and
    // collocation segments generally have dense C matrices
    template <typename Derived>
    struct SegmentDataBase : NumericalBase<Derived>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::SegmentDerived SegmentDerived;
        GALILEO_SEGMENT_TYPEDEF_TEMPLATE(SegmentDerived);

        Scalar cost; //!< Cost integrated over the whole segment

        VectorXs_t C;  //!< Integration constraint value (each row corresponds to the constraint at a node in the segment)
        MatrixXs_t Cx; //!< Jacobian of the integration constraint w.r.t. the state
        MatrixXs_t Cu; //!< Jacobian of the integration constraint w.r.t. the control

    protected:
        // Default constructor: protected
        inline SegmentDataBase()
        {
        }
    };

} // namespace galileo