#ifndef __galileo_predictive_phases_phase_visitors_hxx__
#define __galileo_predictive_phases_phase_visitors_hxx__

#include "galileo/predictive/phases/phase-visitor-base.hpp"
#include "galileo/predictive/phases/phase-visitors.hpp"

namespace galileo
{

    template <typename BasicSpec, typename StateMatrixType, typename ControlParamMatrixType>
    struct PhaseCalcZerothOrderVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCalcZerothOrderVisitor<BasicSpec, StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<const StateMatrixType &, const ControlParamMatrixType &>;

        template <typename PhaseModelType>
        static void algo(
            const PhaseModelBase<PhaseModelType, BasicSpec> &phase_model,
            PhaseDataBase<typename PhaseModelType::Data_t, BasicSpec> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            const Eigen::MatrixBase<ControlParamMatrixType> &ws)
        {
            phase_model.calc(phase_data.derived(), xs.derived(), ws.derived());
        }
    };

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_zeroth_order(
        const PhaseModelTpl<BasicSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws)
    {
        typedef PhaseCalcZerothOrderVisitor<BasicSpec, StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs.derived(), ws.derived()));
    }

    template <typename BasicSpec, typename StateMatrixType, typename ControlParamMatrixType>
    struct PhaseCalcFirstOrderVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCalcFirstOrderVisitor<BasicSpec, StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<const StateMatrixType &, const ControlParamMatrixType &>;

        template <typename PhaseModelType>
        static void algo(
            const PhaseModelBase<PhaseModelType, BasicSpec> &phase_model,
            PhaseDataBase<typename PhaseModelType::Data_t, BasicSpec> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            const Eigen::MatrixBase<ControlParamMatrixType> &ws)
        {
            phase_model.calcDiff(phase_data.derived(), xs.derived(), ws.derived());
        }
    };

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_calc_first_order(
        const PhaseModelTpl<BasicSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        const Eigen::MatrixBase<ControlParamMatrixType> &ws)
    {
        typedef PhaseCalcFirstOrderVisitor<BasicSpec, StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs.derived(), ws.derived()));
    }

    template <typename BasicSpec,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    struct PhaseQuasiStaticVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseQuasiStaticVisitor<BasicSpec, StateMatrixType, ControlParamMatrixType>>
    {
        using ArgsType = boost::fusion::vector<const StateMatrixType &, ControlParamMatrixType &, const int, const typename BasicSpec::NumScalar>;

        template <typename PhaseModelType>
        static void algo(
            const PhaseModelBase<PhaseModelType, BasicSpec> &phase_model,
            PhaseDataBase<typename PhaseModelType::Data_t, BasicSpec> &phase_data,
            const Eigen::MatrixBase<StateMatrixType> &xs,
            Eigen::MatrixBase<ControlParamMatrixType> &ws,
            const int maxiter,
            const typename BasicSpec::NumScalar tol)
        {
            phase_model.quasiStatic(phase_data.derived(), xs.derived(), ws.derived(), maxiter, tol);
        }
    };

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl,
              typename StateMatrixType,
              typename ControlParamMatrixType>
    inline void phase_quasi_static(
        const PhaseModelTpl<BasicSpec, PhaseCollectionTpl> &phase_model,
        PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data,
        const Eigen::MatrixBase<StateMatrixType> &xs,
        Eigen::MatrixBase<ControlParamMatrixType> &ws,
        const int maxiter,
        const typename BasicSpec::NumScalar tol)
    {
        typedef PhaseQuasiStaticVisitor<BasicSpec, StateMatrixType, ControlParamMatrixType> Algo;

        Algo::run(phase_model, phase_data, typename Algo::ArgsType(xs.derived(), ws.derived(), maxiter, tol));
    }

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    struct PhaseCreateDataVisitor
        : fusion::PhaseUnaryVisitorBase<PhaseCreateDataVisitor<BasicSpec, PhaseCollectionTpl>,
                                        PhaseDataTpl<BasicSpec, PhaseCollectionTpl>>
    {
        using PhaseCollection_t = PhaseCollectionTpl<BasicSpec>;
        using PhaseModelVariant_t = PhaseCollection_t::PhaseModelVariant_t;
        using PhaseDataVariant_t = PhaseDataTpl<BasicSpec, PhaseCollectionTpl>;

        template <typename PhaseModelType>
        static PhaseDataVariant_t algo(
            const PhaseModelBase<PhaseModelType, BasicSpec> &phase_model)
        {
            return PhaseDataVariant_t(phase_model.createData());
        }
    };

    template <typename BasicSpec,
              template <typename> class PhaseCollectionTpl>
    inline PhaseDataTpl<BasicSpec, PhaseCollectionTpl> phase_create_data(
        const PhaseModelTpl<BasicSpec, PhaseCollectionTpl> &phase_model)
    {
        typedef PhaseCreateDataVisitor<BasicSpec, PhaseCollectionTpl> Algo;

        return Algo::run(phase_model);
    }

    // Phase data visitors

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseXNextVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNext_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNext_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.XNext_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseXNextVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNext_t phase_XNext_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseXNextVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseXNextxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextx_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.XNextx_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseXNextxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextx_t phase_XNextx_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseXNextxVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseXNextwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextw_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.XNextw_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseXNextwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextw_t phase_XNextw_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseXNextwVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::L_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::L_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.L_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseLVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::L_t phase_L_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseLVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lx_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Lx_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseLxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lx_t phase_Lx_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseLxVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lw_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Lw_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseLwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lw_t phase_Lw_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseLwVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLxxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxx_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Lxx_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseLxxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxx_t phase_Lxx_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseLxxVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLxwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxw_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Lxw_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseLxwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxw_t phase_Lxw_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseLxwVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLwwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lww_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lww_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Lww_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseLwwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lww_t phase_Lww_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseLwwVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseHVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::H_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::H_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.H_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseHVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::H_t phase_H_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseHVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseHxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hx_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Hx_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseHxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hx_t phase_Hx_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseHxVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseHwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hw_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Hw_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseHwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hw_t phase_Hw_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseHwVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseGVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::G_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::G_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.G_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseGVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::G_t phase_G_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseGVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseGxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gx_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Gx_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseGxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gx_t phase_Gx_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseGxVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseGwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gw_t;
        using ArgsType = boost::fusion::vector<const int>;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data, const int i) const
        {
            return phase_data.Gw_at_seg_i(i);
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
        {
            return boost::apply_visitor(PhaseGwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data, i);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gw_t phase_Gw_at_seg_i(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data, const int i)
    {
        using Algo = PhaseGwVisitor<BasicSpec, PhaseCollectionTpl>;
        return Algo::run(phase_data, typename Algo::ArgsType(i));
    }

} // namespace galileo

#endif // __galileo_predictive_phases_phase_visitors_hxx__
