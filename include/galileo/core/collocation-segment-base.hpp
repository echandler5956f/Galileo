#pragma once

#include "galileo/core/segment-base.hpp"

#define GALILEO_SHOOTING_SEGMENT_TYPEDEF_GENERIC(CollocationSegment, TYPENAME)                                  \
    typedef TYPENAME traits<CollocationSegment>::Scalar Scalar;                                                 \
    typedef TYPENAME traits<CollocationSegment>::CollocationSegmentModelDerived CollocationSegmentModelDerived; \
    typedef TYPENAME traits<CollocationSegment>::CollocationSegmentDataDerived CollocationSegmentDataDerived;   \
    typedef TYPENAME traits<CollocationSegment>::SegmentModelDerived SegmentModelDerived;                       \
    typedef TYPENAME traits<CollocationSegment>::SegmentDataDerived SegmentDataDerived;                         \
    typedef TYPENAME traits<CollocationSegment>::NodeModelDerived NodeModelDerived;                             \
    typedef TYPENAME traits<CollocationSegment>::NodeDataDerived NodeDataDerived;                               \
    typedef TYPENAME traits<CollocationSegment>::StateModel StateModel;                                         \
    typedef TYPENAME traits<CollocationSegment>::ControlModel ControlModel;                                     \
    typedef TYPENAME traits<CollocationSegment>::MathBase MathBase;                                             \
    typedef TYPENAME MathBase::VectorXs VectorXs;                                                               \
    typedef TYPENAME MathBase::MatrixXs MatrixXs;

#define GALILEO_SHOOTING_SEGMENT_TYPEDEF_TEMPLATE(CollocationSegment) \
    GALILEO_SHOOTING_SEGMENT_TYPEDEF_GENERIC(CollocationSegment, typename)

#define GALILEO_SHOOTING_SEGMENT_CAST_TYPE_SPECIALIZATION(CollocationSegmentModelTpl) \
    template <typename Scalar, typename NewScalar>                                    \
    struct CastType<NewScalar, CollocationSegmentModelTpl<Scalar>>                    \
    {                                                                                 \
        typedef CollocationSegmentModelTpl<NewScalar> type;                           \
    }

namespace galileo
{
    template <typename Derived>
    class CollocationSegmentModelBase : SegmentModelBase<CollocationSegmentModelBase<Derived>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::CollocationSegmentDerived CollocationSegmentDerived;
        GALILEO_SHOOTING_SEGMENT_TYPEDEF_TEMPLATE(CollocationSegmentDerived);

        CollocationSegmentModelDerived &derived()
        {
            return *static_cast<Derived *>(this);
        }

        const CollocationSegmentModelDerived &derived() const
        {
            return *static_cast<const Derived *>(this);
        }

        void calc(const boost::shared_ptr<CollocationSegmentDataDerived> &data,
                  const Eigen::Ref<const MatrixXs> &x,
                  const Eigen::Ref<const MatrixXs> &u)
        {
            derived().calc(data, x, u);
        }

        void calc(const boost::shared_ptr<CollocationSegmentDataDerived> &data,
                  const Eigen::Ref<const MatrixXs> &x)
        {
            derived().calc(data, x);
        }

        void calcDiff(const boost::shared_ptr<CollocationSegmentDataDerived> &data,
                      const Eigen::Ref<const MatrixXs> &x,
                      const Eigen::Ref<const MatrixXs> &u)
        {
            derived().calcDiff(data, x, u);
        }

        void calcDiff(const boost::shared_ptr<CollocationSegmentDataDerived> &data,
                      const Eigen::Ref<const MatrixXs> &x)
        {
            derived().calcDiff(data, x);
        }

        void quasiStatic(const boost::shared_ptr<CollocationSegmentDataDerived> &data,
                         Eigen::Ref<MatrixXs> u, const Eigen::Ref<const MatrixXs> &x,
                         const std::size_t maxiter = 100, const Scalar tol = Scalar(1e-9))
        {
            derived().quasiStatic(data, u, x, maxiter, tol);
        }

        std::vector<VectorXs> quasiStatic_xs(const boost::shared_ptr<CollocationSegmentDataDerived> &data,
                                             const Eigen::Ref<const MatrixXs> &x,
                                             const std::size_t maxiter = 100,
                                             const Scalar tol = Scalar(1e-9))
        {
            return derived().quasiStatic_xs(data, x, maxiter, tol);
        }

        void circularAppend(boost::shared_ptr<NodeModelDerived> model, boost::shared_ptr<NodeDataDerived> data)
        {
            derived().circularAppend(model, data);
        }

        void circularAppend(boost::shared_ptr<NodeModelDerived> model)
        {
            derived().circularAppend(model);
        }

        void updateNode(const std::size_t &id, boost::shared_ptr<NodeModelDerived> model, boost::shared_ptr<NodeDataDerived> data)
        {
            derived().updateNode(id, model, data);
        }

        void updateNode(const std::size_t &id, boost::shared_ptr<NodeModelDerived> model)
        {
            derived().updateNode(id, model);
        }

        std::size_t getNumNodes() const
        {
            return derived().getNumNodes();
        }

        const std::vector<boost::shared_ptr<NodeModelDerived>> &getNodes() const
        {
            return derived().getNodes();
        }

        const std::vector<boost::shared_ptr<NodeDataDerived>> &getNodesData() const
        {
            return derived().getNodesData();
        }

        const ControlModel &getControl() const
        {
            return derived().getControl();
        }

        const VectorXs &getNodeTimes() const
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
        // Prevent the construction of stand-alone CollocationSegmentModelBase.
        inline CollocationSegmentModelBase() : Ns_(0), h_(0), nx_(0), ndx_(0), nu_(0), is_updated_(false)
        {
        }

        // Copy constructor: protected
        // Copy of stand-alone CollocationSegmentModelBase are prevented, but can be used from inheriting
        // objects.
        inline CollocationSegmentModelBase(const CollocationSegmentModelBase &clone)
        {
            *this = clone;
        }

        // Copy operator: protected
        // Copy of stand-alone CollocationSegmentModelBase are prevented, but can be used from inheriting
        // objects.
        inline CollocationSegmentModelBase &operator=(const CollocationSegmentModelBase &clone)
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

        using SegmentModelBase<CollocationSegmentModel<Derived>>::Ns_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::nodes_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::nodes_data_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::control_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::node_times_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::h_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::nx_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::ndx_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::nu_;
        using SegmentModelBase<CollocationSegmentModel<Derived>>::is_updated_;
    };

    template <typename Derived>
    struct CollocationSegmentData : SegmentDataBase<CollocationSegmentData<Derived>>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::CollocationSegmentDerived CollocationSegmentDerived;
        GALILEO_SHOOTING_SEGMENT_TYPEDEF_TEMPLATE(CollocationSegmentDerived);

        using SegmentDataBase<CollocationSegmentData<Derived>>::cost;

        using SegmentDataBase<CollocationSegmentData<Derived>>::C;
        using SegmentDataBase<CollocationSegmentData<Derived>>::Cx;
        using SegmentDataBase<CollocationSegmentData<Derived>>::Cu;

    protected:
        // Default constructor: protected
        inline CollocationSegmentData()
        {
        }
    };
}