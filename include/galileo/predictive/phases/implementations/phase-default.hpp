#ifndef __galileo_predictive_phases_phase_default_hpp__
#define __galileo_predictive_phases_phase_default_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename PhaseSpec>
    struct PhaseDefaultTpl;

    template <typename PhaseSpec>
    struct traits<PhaseDefaultTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        // GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseDefaultTpl<PS>;
        using Model_t = PhaseModelDefaultTpl<PS>;
        using Data_t = PhaseDataDefaultTpl<PS>;
    };

    template <typename PhaseSpec>
    struct traits<PhaseDataDefaultTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseDefaultTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct traits<PhaseModelDefaultTpl<PhaseSpec>>
    {
        using PS = PhaseSpec;

        using Meta_t = PhaseDefaultTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec>
    struct PhaseDataDefaultTpl
        : public PhaseDataBase<PhaseDataDefaultTpl<PhaseSpec>, typename PhaseSpec::BS>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseDefaultTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseDataBase<PhaseDataDefaultTpl<PS>, typename PS::BS>;

        PhaseDataDefaultTpl(const Model_t &model)
            : Base()
        {
            segments.reserve(model.get_segments().size());
            for (const auto &segment : model.get_segments())
            {
                segments.push_back(segment.createData());
            }
        }

        std::vector<SegmentData_t> &get_segments()
        {
            return segments;
        }
        const std::vector<SegmentData_t> &get_segments() const
        {
            return segments;
        }

        XNext_t &XNext_at_i_accessor(const int i)
        {
            return segments[i].XNext();
        }
        const XNext_t &XNext_at_i_accessor(const int i) const
        {
            return segments[i].XNext();
        }

        XNextx_t &XNextx_at_i_accessor(const int i)
        {
            return segments[i].XNextx();
        }
        const XNextx_t &XNextx_at_i_accessor(const int i) const
        {
            return segments[i].XNextx();
        }

        XNextw_t &XNextw_at_i_accessor(const int i)
        {
            return segments[i].XNextw();
        }
        const XNextw_t &XNextw_at_i_accessor(const int i) const
        {
            return segments[i].XNextw();
        }

        L_t &L_at_i_accessor(const int i)
        {
            return segments[i].L();
        }
        const L_t &L_at_i_accessor(const int i) const
        {
            return segments[i].L();
        }

        Lx_t &Lx_at_i_accessor(const int i)
        {
            return segments[i].Lx();
        }
        const Lx_t &Lx_at_i_accessor(const int i) const
        {
            return segments[i].Lx();
        }

        Lw_t &Lw_at_i_accessor(const int i)
        {
            return segments[i].Lw();
        }
        const Lw_t &Lw_at_i_accessor(const int i) const
        {
            return segments[i].Lw();
        }

        Lxx_t &Lxx_at_i_accessor(const int i)
        {
            return segments[i].Lxx();
        }
        const Lxx_t &Lxx_at_i_accessor(const int i) const
        {
            return segments[i].Lxx();
        }

        Lxw_t &Lxw_at_i_accessor(const int i)
        {
            return segments[i].Lxw();
        }
        const Lxw_t &Lxw_at_i_accessor(const int i) const
        {
            return segments[i].Lxw();
        }

        Lww_t &Lww_at_i_accessor(const int i)
        {
            return segments[i].Lww();
        }
        const Lww_t &Lww_at_i_accessor(const int i) const
        {
            return segments[i].Lww();
        }

        H_t &H_at_i_accessor(const int i)
        {
            return segments[i].H();
        }
        const H_t &H_at_i_accessor(const int i) const
        {
            return segments[i].H();
        }

        Hx_t &Hx_at_i_accessor(const int i)
        {
            return segments[i].Hx();
        }
        const Hx_t &Hx_at_i_accessor(const int i) const
        {
            return segments[i].Hx();
        }

        Hw_t &Hw_at_i_accessor(const int i)
        {
            return segments[i].Hw();
        }
        const Hw_t &Hw_at_i_accessor(const int i) const
        {
            return segments[i].Hw();
        }

        G_t &G_at_i_accessor(const int i)
        {
            return segments[i].G();
        }
        const G_t &G_at_i_accessor(const int i) const
        {
            return segments[i].G();
        }

        Gx_t &Gx_at_i_accessor(const int i)
        {
            return segments[i].Gx();
        }
        const Gx_t &Gx_at_i_accessor(const int i) const
        {
            return segments[i].Gx();
        }

        Gw_t &Gw_at_i_accessor(const int i)
        {
            return segments[i].Gw();
        }
        const Gw_t &Gw_at_i_accessor(const int i) const
        {
            return segments[i].Gw();
        }

        std::vector<SegmentData_t> segments;

    }; // struct PhaseDataDefaultTpl

    template <typename PhaseSpec>
    class PhaseModelDefaultTpl
        : public PhaseModelBase<PhaseModelDefaultTpl<PhaseSpec>, typename PhaseSpec::BS>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseDefaultTpl<PS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseModelBase<PhaseModelDefaultTpl<PS>, typename PS::BS>;

        PhaseModelDefaultTpl(const PS &ps)
            : Base(), ps_(ps)
        {
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs,
                  const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            int i = 0;
            for (auto model_it = segments_.begin(), data_it = data.segments.begin();
                 model_it != segments_.end(); ++model_it, ++data_it, ++i)
            {
                model_it->calc(*data_it, col(xs, i), col(ws, i));
            }
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calcDiff(Data_t &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            int i = 0;
            for (auto model_it = segments_.begin(), data_it = data.segments.begin();
                 model_it != segments_.end(); ++model_it, ++data_it, ++i)
            {
                model_it->calcDiff(*data_it, col(xs, i), col(ws, i));
            }
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void quasiStatic(Data_t &data, const Eigen::MatrixBase<StateMatrixType> &xs,
                         Eigen::MatrixBase<ControlParamMatrixType> &ws,
                         const int maxiter, const NumScalar &tol) const
        {
            int i = 0;
            for (auto model_it = segments_.begin(), data_it = data.segments.begin();
                 model_it != segments_.end(); ++model_it, ++data_it, ++i)
            {
                model_it->quasiStatic(*data_it, col(xs, i), col(ws, i), maxiter, tol);
            }
        }

        Data_t createData() const
        {
            return Data_t(*this);
        }

        const PS &get_ps() const
        {
            return ps_.get();
        }

        const std::vector<SegmentModel_t> &get_segments() const
        {
            return segments_;
        }

    protected:
        std::reference_wrapper<const PS> ps_;
        std::vector<SegmentModel_t> segments_;

    }; // class PhaseModelDefaultTpl

} // namespace galileo

#endif // __galileo_predictive_phases_phase_default_hpp__
