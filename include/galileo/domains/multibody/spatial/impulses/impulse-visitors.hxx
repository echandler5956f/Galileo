#ifndef __galileo_multibody_spatial_impulses_impulse_visitors_hxx__
#define __galileo_multibody_spatial_impulses_impulse_visitors_hxx__

#include "galileo/domains/multibody/spatial/impulses/impulse-visitor-base.hpp"
#include "galileo/domains/multibody/spatial/impulses/impulse-visitors.hpp"

namespace galileo
{

    // Impulse model visitors

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename StateVectorType>
    struct ImpulseCalcZerothOrderVisitor
        : fusion::ImpulseUnaryVisitorBase<
              ImpulseCalcZerothOrderVisitor<PhaseSpec, ImpulseCollectionTpl, StateVectorType>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &>;

        template <typename ImpulseModelType>
        static void algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                         ImpulseDataBase<typename ImpulseModelType::Data_t, PhaseSpec> &impulse_data,
                         const Eigen::MatrixBase<StateVectorType> &x)
        {
            impulse_model.calc(impulse_data.derived(), x.derived());
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename StateVectorType>
    inline void impulse_calc_zeroth_order(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                          ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                                          const Eigen::MatrixBase<StateVectorType> &x)
    {
        typedef ImpulseCalcZerothOrderVisitor<PhaseSpec, ImpulseCollectionTpl, StateVectorType> Algo;

        Algo::run(impulse_model, impulse_data, typename Algo::ArgsType(x.derived()));
    }

    template <typename PhaseSpec, template <typename PS> class ImpulseCollectionTpl, typename StateVectorType>
    struct ImpulseCalcFirstOrderVisitor
        : fusion::ImpulseUnaryVisitorBase<
              ImpulseCalcFirstOrderVisitor<PhaseSpec, ImpulseCollectionTpl, StateVectorType>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &>;

        template <typename ImpulseModelType>
        static void algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                         ImpulseDataBase<typename ImpulseModelType::Data_t, PhaseSpec> &impulse_data,
                         const Eigen::MatrixBase<StateVectorType> &x)
        {
            impulse_model.calcDiff(impulse_data.derived(), x.derived());
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename StateVectorType>
    inline void impulse_calc_first_order(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                         ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                                         const Eigen::MatrixBase<StateVectorType> &x)
    {
        typedef ImpulseCalcFirstOrderVisitor<PhaseSpec, ImpulseCollectionTpl, StateVectorType> Algo;

        Algo::run(impulse_model, impulse_data, typename Algo::ArgsType(x.derived()));
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename ForceVectorType>
    struct ImpulseUpdateForceVisitor
        : fusion::ImpulseUnaryVisitorBase<ImpulseUpdateForceVisitor<PhaseSpec, ImpulseCollectionTpl, ForceVectorType>>
    {
        using ArgsType = boost::fusion::vector<const ForceVectorType &>;

        template <typename ImpulseModelType>
        static void algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                         ImpulseDataBase<typename ImpulseModelType::Data_t, PhaseSpec> &impulse_data,
                         const Eigen::MatrixBase<ForceVectorType> &force)
        {
            impulse_model.updateForce(impulse_data.derived(), force.derived());
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename ForceVectorType>
    inline void impulse_update_force(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                     ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                                     const Eigen::MatrixBase<ForceVectorType> &force)
    {
        typedef ImpulseUpdateForceVisitor<PhaseSpec, ImpulseCollectionTpl, ForceVectorType> Algo;

        Algo::run(impulse_model, impulse_data, typename Algo::ArgsType(force.derived()));
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename MatrixNcNdxType>
    struct ImpulseUpdateForceDiffVisitor
        : fusion::ImpulseUnaryVisitorBase<
              ImpulseUpdateForceDiffVisitor<PhaseSpec, ImpulseCollectionTpl, MatrixNcNdxType>>
    {
        using ArgsType = boost::fusion::vector<const MatrixNcNdxType &>;

        template <typename ImpulseModelType>
        static void algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                         ImpulseDataBase<typename ImpulseModelType::Data_t, PhaseSpec> &impulse_data,
                         const Eigen::MatrixBase<MatrixNcNdxType> &df_dx)
        {
            impulse_model.updateForceDiff(impulse_data.derived(), df_dx.derived());
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl, typename MatrixNcNdxType>
    inline void impulse_update_force_diff(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                          ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data,
                                          const Eigen::MatrixBase<MatrixNcNdxType> &df_dx)
    {
        typedef ImpulseUpdateForceDiffVisitor<PhaseSpec, ImpulseCollectionTpl, MatrixNcNdxType> Algo;

        Algo::run(impulse_model, impulse_data, typename Algo::ArgsType(df_dx.derived()));
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseSetZeroForceVisitor
        : fusion::ImpulseUnaryVisitorBase<ImpulseSetZeroForceVisitor<PhaseSpec, ImpulseCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<>;

        template <typename ImpulseModelType>
        static void algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                         ImpulseDataBase<typename ImpulseModelType::Data_t, PhaseSpec> &impulse_data)
        {
            impulse_model.setZeroForce(impulse_data.derived());
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline void impulse_set_zero_force(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                       ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        typedef ImpulseSetZeroForceVisitor<PhaseSpec, ImpulseCollectionTpl> Algo;

        Algo::run(impulse_model, impulse_data, typename Algo::ArgsType());
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseSetZeroForceDiffVisitor
        : fusion::ImpulseUnaryVisitorBase<ImpulseSetZeroForceDiffVisitor<PhaseSpec, ImpulseCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<>;

        template <typename ImpulseModelType>
        static void algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                         ImpulseDataBase<typename ImpulseModelType::Data_t, PhaseSpec> &impulse_data)
        {
            impulse_model.setZeroForceDiff(impulse_data.derived());
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline void impulse_set_zero_force_diff(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                                            ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        typedef ImpulseSetZeroForceDiffVisitor<PhaseSpec, ImpulseCollectionTpl> Algo;

        Algo::run(impulse_model, impulse_data, typename Algo::ArgsType());
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseCreateDataVisitor
        : fusion::ImpulseUnaryVisitorBase<ImpulseCreateDataVisitor<PhaseSpec, ImpulseCollectionTpl>,
                                          ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>>
    {
        using ArgsType = boost::fusion::vector<typename PhaseSpec::RobotData_t *const>;
        using ImpulseCollection_t = ImpulseCollectionTpl<PhaseSpec>;
        using ImpulseModelVariant_t = ImpulseCollection_t::ImpulseModelVariant_t;
        using ImpulseDataVariant_t = ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>;

        template <typename ImpulseModelType>
        static ImpulseDataVariant_t algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                                         typename PhaseSpec::RobotData_t *const robot)
        {
            return ImpulseDataVariant_t(impulse_model.createData(robot));
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> impulse_create_data(
        const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
        typename PhaseSpec::RobotData_t *const robot)
    {
        typedef ImpulseCreateDataVisitor<PhaseSpec, ImpulseCollectionTpl> Algo;
        return Algo::run(impulse_model, typename Algo::ArgsType(robot));
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseGetIdVisitor
        : boost::static_visitor<typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t>
    {
        using ReturnType = typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t;

        template <typename ImpulseModelType>
        ReturnType operator()(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model) const
        {
            return impulse_model.get_id();
        }

        static ReturnType run(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model)
        {
            return boost::apply_visitor(ImpulseGetIdVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_model);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t impulse_get_id(
        const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model)
    {
        return ImpulseGetIdVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_model);
    }

    template <typename PhaseSpec, template <typename PS> class ImpulseCollectionTpl>
    struct ImpulseSetIdVisitor : boost::static_visitor<void>
    {
        using ArgsType =
            boost::fusion::vector<const typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t &>;

        template <typename ImpulseModelType>
        static void algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                         const typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t &id)
        {
            impulse_model.set_id(id);
        }
    };

    template <typename PhaseSpec, template <typename PS> class ImpulseCollectionTpl>
    inline void impulse_set_id(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
                               const typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t &id)
    {
        typedef ImpulseSetIdVisitor<PhaseSpec, ImpulseCollectionTpl> Algo;

        Algo::run(impulse_model, typename Algo::ArgsType(id));
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseGetTypeVisitor
        : boost::static_visitor<typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t>
    {
        using ReturnType = typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t;

        template <typename ImpulseModelType>
        ReturnType operator()(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model) const
        {
            return impulse_model.get_type();
        }

        static ReturnType run(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model)
        {
            return boost::apply_visitor(ImpulseGetTypeVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_model);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t impulse_get_type(
        const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model)
    {
        return ImpulseGetTypeVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_model);
    }

    template <typename PhaseSpec, template <typename PS> class ImpulseCollectionTpl>
    struct ImpulseSetTypeVisitor : boost::static_visitor<void>
    {
        using ArgsType =
            boost::fusion::vector<const typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t &>;

        template <typename ImpulseModelType>
        static void algo(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model,
                         const typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t &type)
        {
            impulse_model.set_type(type);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline void impulse_set_type(
        const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model,
        const typename ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t &type)
    {
        typedef ImpulseSetTypeVisitor<PhaseSpec, ImpulseCollectionTpl> Algo;

        Algo::run(impulse_model, typename Algo::ArgsType(type));
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseGetNcVisitor : boost::static_visitor<int>
    {
        template <typename ImpulseModelType>
        int operator()(const ImpulseModelBase<ImpulseModelType, PhaseSpec> &impulse_model) const
        {
            return impulse_model.get_nc();
        }

        static int run(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model)
        {
            return boost::apply_visitor(ImpulseGetNcVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_model);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline int impulse_get_nc(const ImpulseModelTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_model)
    {
        return ImpulseGetNcVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_model);
    }

    // Impulse data visitors

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseRobotDataVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::RobotData_t>
    {

        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::RobotDataPointer_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.robot();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseRobotDataVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::RobotData_t *impulse_robot_data(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseRobotDataVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseFrameVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.frame();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseFrameVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::FrameIndex_t impulse_frame(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseFrameVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseTypeDataVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.type();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseTypeDataVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::ReferenceFrame_t impulse_type_data(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseTypeDataVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseJmFVisitor : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::SE3_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::SE3_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.jMf();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseJmFVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::SE3_t impulse_jMf(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseJmFVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseJcVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNv_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNv_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.Jc();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseJcVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNv_t impulse_Jc(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseJcVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseFVisitor : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::Force_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::Force_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.f();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseFVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::Force_t impulse_f(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseFVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseFextVisitor : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::Force_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::Force_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.fext();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseFextVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::Force_t impulse_fext(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseFextVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulsedFdXVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNdx_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNdx_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.df_dx();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulsedFdXVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNdx_t impulse_df_dx(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulsedFdXVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulsedFdUVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNu_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNu_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.df_du();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulsedFdUVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNu_t impulse_df_du(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulsedFdUVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseFXjVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::ActionMatrix_t>
    {

        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::ActionMatrix_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.fXj();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseFXjVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::ActionMatrix_t impulse_fXj(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseFXjVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseDV0dQVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNv_t>
    {

        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNv_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.dv0_dq();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseDV0dQVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNcNv_t impulse_dv0_dq(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseDV0dQVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    struct ImpulseDtauDqVisitor
        : boost::static_visitor<typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNv_t>
    {
        using ReturnType = typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNv_t;

        template <typename ImpulseDataType>
        ReturnType operator()(const ImpulseDataBase<ImpulseDataType, PhaseSpec> &impulse_data) const
        {
            return impulse_data.dtau_dq();
        }

        static ReturnType run(const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
        {
            return boost::apply_visitor(ImpulseDtauDqVisitor<PhaseSpec, ImpulseCollectionTpl>(), impulse_data);
        }
    };

    template <typename PhaseSpec, template <typename> class ImpulseCollectionTpl>
    inline typename ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl>::MatrixNv_t impulse_dtau_dq(
        const ImpulseDataTpl<PhaseSpec, ImpulseCollectionTpl> &impulse_data)
    {
        return ImpulseDtauDqVisitor<PhaseSpec, ImpulseCollectionTpl>::run(impulse_data);
    }

} // namespace galileo

#endif // __galileo_multibody_spatial_impulses_impulse_visitors_hxx__
