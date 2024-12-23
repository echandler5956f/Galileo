#pragma once

#include "galileo/core/segment-base.hpp"

#define GALILEO_SHOOTING_SEGMENT_TYPEDEF_GENERIC(ShootingSegment, TYPENAME)                            \
    typedef TYPENAME traits<ShootingSegment>::Scalar Scalar;                                           \
    typedef TYPENAME traits<ShootingSegment>::ShootingSegmentModelDerived ShootingSegmentModelDerived; \
    typedef TYPENAME traits<ShootingSegment>::ShootingSegmentDataDerived ShootingSegmentDataDerived;   \
    typedef TYPENAME traits<ShootingSegment>::SegmentModelDerived SegmentModelDerived;                 \
    typedef TYPENAME traits<ShootingSegment>::SegmentDataDerived SegmentDataDerived;                   \
    typedef TYPENAME traits<ShootingSegment>::NodeModelDerived NodeModelDerived;                       \
    typedef TYPENAME traits<ShootingSegment>::NodeDataDerived NodeDataDerived;                         \
    typedef TYPENAME traits<ShootingSegment>::StateModel StateModel;                                   \
    typedef TYPENAME traits<ShootingSegment>::ControlModel ControlModel;                               \
    typedef TYPENAME traits<ShootingSegment>::MathBase MathBase;                                       \
    typedef TYPENAME MathBase::VectorXs VectorXs;                                                      \
    typedef TYPENAME MathBase::MatrixXs MatrixXs;

#define GALILEO_SHOOTING_SEGMENT_TYPEDEF_TEMPLATE(ShootingSegment) \
    GALILEO_SHOOTING_SEGMENT_TYPEDEF_GENERIC(ShootingSegment, typename)

#define GALILEO_SHOOTING_SEGMENT_CAST_TYPE_SPECIALIZATION(ShootingSegmentModelTpl) \
    template <typename Scalar, typename NewScalar>                                 \
    struct CastType<NewScalar, ShootingSegmentModelTpl<Scalar>>                    \
    {                                                                              \
        typedef ShootingSegmentModelTpl<NewScalar> type;                           \
    }

namespace galileo
{
    template <typename Derived>
    class ShootingSegmentModelBase : SegmentModelBase<ShootingSegmentModelBase<Derived>>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::ShootingSegmentDerived ShootingSegmentDerived;
        GALILEO_SHOOTING_SEGMENT_TYPEDEF_TEMPLATE(ShootingSegmentDerived);

        ShootingSegmentModelDerived &derived()
        {
            return *static_cast<Derived *>(this);
        }

        const ShootingSegmentModelDerived &derived() const
        {
            return *static_cast<const Derived *>(this);
        }

        void calc(const boost::shared_ptr<ShootingSegmentDataDerived> &data,
                  const Eigen::Ref<const MatrixXs> &x,
                  const Eigen::Ref<const MatrixXs> &u)
        {
            derived().calc(data, x, u);
        }

        void calc(const boost::shared_ptr<ShootingSegmentDataDerived> &data,
                  const Eigen::Ref<const MatrixXs> &x)
        {
            derived().calc(data, x);
        }

        void calcDiff(const boost::shared_ptr<ShootingSegmentDataDerived> &data,
                      const Eigen::Ref<const MatrixXs> &x,
                      const Eigen::Ref<const MatrixXs> &u)
        {
            derived().calcDiff(data, x, u);
        }

        void calcDiff(const boost::shared_ptr<ShootingSegmentDataDerived> &data,
                      const Eigen::Ref<const MatrixXs> &x)
        {
            derived().calcDiff(data, x);
        }

        void quasiStatic(const boost::shared_ptr<ShootingSegmentDataDerived> &data,
                         Eigen::Ref<MatrixXs> u, const Eigen::Ref<const MatrixXs> &x,
                         const std::size_t maxiter = 100, const Scalar tol = Scalar(1e-9))
        {
            derived().quasiStatic(data, u, x, maxiter, tol);
        }

        std::vector<VectorXs> quasiStatic_xs(const boost::shared_ptr<ShootingSegmentDataDerived> &data,
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
        // Prevent the construction of stand-alone ShootingSegmentModelBase.
        inline ShootingSegmentModelBase() : Ns_(0), h_(0), nx_(0), ndx_(0), nu_(0), is_updated_(false)
        {
        }

        // Copy constructor: protected
        // Copy of stand-alone ShootingSegmentModelBase are prevented, but can be used from inheriting
        // objects.
        inline ShootingSegmentModelBase(const ShootingSegmentModelBase &clone)
        {
            *this = clone;
        }

        // Copy operator: protected
        // Copy of stand-alone ShootingSegmentModelBase are prevented, but can be used from inheriting
        // objects.
        inline ShootingSegmentModelBase &operator=(const ShootingSegmentModelBase &clone)
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

        using SegmentModelBase<ShootingSegmentModel<Derived>>::Ns_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::nodes_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::nodes_data_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::control_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::node_times_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::h_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::nx_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::ndx_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::nu_;
        using SegmentModelBase<ShootingSegmentModel<Derived>>::is_updated_;
    };

    template <typename Derived>
    struct ShootingSegmentData : SegmentDataBase<ShootingSegmentData<Derived>>
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef typename traits<Derived>::ShootingSegmentDerived ShootingSegmentDerived;
        GALILEO_SHOOTING_SEGMENT_TYPEDEF_TEMPLATE(ShootingSegmentDerived);

        using SegmentDataBase<ShootingSegmentData<Derived>>::cost;

        using SegmentDataBase<ShootingSegmentData<Derived>>::C;
        using SegmentDataBase<ShootingSegmentData<Derived>>::Cx;
        using SegmentDataBase<ShootingSegmentData<Derived>>::Cu;

    protected:
        // Default constructor: protected
        inline ShootingSegmentData()
        {
        }
    };
}