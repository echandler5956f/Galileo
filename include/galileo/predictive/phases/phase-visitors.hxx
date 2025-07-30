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

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.XNext();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseXNextVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNext_t phase_XNext(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseXNextVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseXNextxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextx_t;

        template <typename PhaseDataDerived>
        ReturnType operator()(const PhaseDataBase<PhaseDataDerived, BasicSpec> &phase_data) const
        {
            return phase_data.XNextx();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseXNextxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextx_t phase_XNextx(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseXNextxVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseXNextwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextw_t;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.XNextw();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseXNextwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::XNextw_t phase_XNextw(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseXNextwVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::L_t>
    {

        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::L_t;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.L();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseLVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::L_t phase_L(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseLVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lx_t;

        template <typename PhaseDataDerived>
        ReturnType operator()(const PhaseDataBase<PhaseDataDerived, BasicSpec> &phase_data) const
        {
            return phase_data.Lx();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseLxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lx_t phase_Lx(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseLxVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lw_t;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.Lw();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseLwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lw_t phase_Lw(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseLwVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLxxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxx_t;

        template <typename PhaseDataDerived>
        ReturnType operator()(const PhaseDataBase<PhaseDataDerived, BasicSpec> &phase_data) const
        {
            return phase_data.Lxx();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseLxxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxx_t phase_Lxx(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseLxxVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLxwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxw_t;

        template <typename PhaseDataDerived>
        ReturnType operator()(const PhaseDataBase<PhaseDataDerived, BasicSpec> &phase_data) const
        {
            return phase_data.Lxw();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseLxwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lxw_t phase_Lxw(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseLxwVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseLwwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lww_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lww_t;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.Lww();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseLwwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Lww_t phase_Lww(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseLwwVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseHVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::H_t>
    {

        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::H_t;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.H();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseHVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::H_t phase_H(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseHVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseHxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hx_t;

        template <typename PhaseDataDerived>
        ReturnType operator()(const PhaseDataBase<PhaseDataDerived, BasicSpec> &phase_data) const
        {
            return phase_data.Hx();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseHxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hx_t phase_Hx(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseHxVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseHwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hw_t;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.Hw();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseHwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Hw_t phase_Hw(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseHwVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseGVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::G_t>
    {

        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::G_t;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.G();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseGVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::G_t phase_G(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseGVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseGxVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gx_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gx_t;

        template <typename PhaseDataDerived>
        ReturnType operator()(const PhaseDataBase<PhaseDataDerived, BasicSpec> &phase_data) const
        {
            return phase_data.Gx();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseGxVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gx_t phase_Gx(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseGxVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    struct PhaseGwVisitor
        : boost::static_visitor<typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gw_t>
    {
        using ReturnType = typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gw_t;

        template <typename PhaseDataType>
        ReturnType operator()(const PhaseDataBase<PhaseDataType, BasicSpec> &phase_data) const
        {
            return phase_data.Gw();
        }

        static ReturnType run(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
        {
            return boost::apply_visitor(PhaseGwVisitor<BasicSpec, PhaseCollectionTpl>(), phase_data);
        }
    };

    template <typename BasicSpec, template <typename> class PhaseCollectionTpl>
    inline typename PhaseDataTpl<BasicSpec, PhaseCollectionTpl>::Gw_t phase_Gw(const PhaseDataTpl<BasicSpec, PhaseCollectionTpl> &phase_data)
    {
        return PhaseGwVisitor<BasicSpec, PhaseCollectionTpl>::run(phase_data);
    }

} // namespace galileo

#endif // __galileo_predictive_phases_phase_visitors_hxx__
