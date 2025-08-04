#ifndef __galileo_predictive_phases_phase_default_hpp__
#define __galileo_predictive_phases_phase_default_hpp__

#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-spec.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class JumpTpl>
    struct PhaseDefaultTpl;

    template <typename PhaseSpec,
              template <typename PS> class JumpTpl>
    struct traits<PhaseDefaultTpl<PhaseSpec, JumpTpl>>
    {
        using PS = PhaseSpec;
        using SpecOfBaseClass = typename PS::BS;

        using Meta_t = PhaseDefaultTpl<PS, JumpTpl>;
        using Model_t = PhaseModelDefaultTpl<PS, JumpTpl>;
        using Data_t = PhaseDataDefaultTpl<PS, JumpTpl>;

        using JumpMeta_t = JumpTpl<PS>;
        using JumpModel_t = typename traits<JumpMeta_t>::Model_t;
        using JumpData_t = typename traits<JumpMeta_t>::Data_t;

        // Add the missing types from PhaseSpec
        using XNext_t = typename PS::XNext_t;
        using XNextx_t = typename PS::XNextx_t;
        using XNextw_t = typename PS::XNextw_t;
        using L_t = typename PS::L_t;
        using Lx_t = typename PS::Lx_t;
        using Lw_t = typename PS::Lw_t;
        using Lxx_t = typename PS::Lxx_t;
        using Lxw_t = typename PS::Lxw_t;
        using Lww_t = typename PS::Lww_t;
        using H_t = typename PS::H_t;
        using Hx_t = typename PS::Hx_t;
        using Hw_t = typename PS::Hw_t;
        using G_t = typename PS::G_t;
        using Gx_t = typename PS::Gx_t;
        using Gw_t = typename PS::Gw_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class JumpTpl>
    struct traits<PhaseDataDefaultTpl<PhaseSpec, JumpTpl>>
    {
        using SpecOfBaseClass = typename PhaseSpec::BS;
        using Meta_t = PhaseDefaultTpl<PhaseSpec, JumpTpl>;
    };

    template <typename PhaseSpec,
              template <typename PS> class JumpTpl>
    struct traits<PhaseModelDefaultTpl<PhaseSpec, JumpTpl>>
    {
        using SpecOfBaseClass = typename PhaseSpec::BS;
        using Meta_t = PhaseDefaultTpl<PhaseSpec, JumpTpl>;
    };

    template <typename PhaseSpec,
              template <typename PS> class JumpTpl>
    struct PhaseDataDefaultTpl
        : public PhaseDataBase<PhaseDataDefaultTpl<PhaseSpec, JumpTpl>, typename PhaseSpec::BS>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseDefaultTpl<PS, JumpTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseDataBase<PhaseDataDefaultTpl<PS, JumpTpl>, typename PS::BS>;

        using JumpMeta_t = typename traits<Meta_t>::JumpMeta_t;
        using JumpModel_t = typename traits<Meta_t>::JumpModel_t;
        using JumpData_t = typename traits<Meta_t>::JumpData_t;

        PhaseDataDefaultTpl(const Model_t &model)
            : Base(),
              jump(model.get_jump().createData())
        {
            segments.reserve(model.get_segments().size());
            for (const auto &segment : model.get_segments())
            {
                segments.push_back(segment.createData());
            }
        }

        SegmentDataVector_t &get_segments()
        {
            return segments;
        }
        const SegmentDataVector_t &get_segments() const
        {
            return segments;
        }

        XNext_t &XNext_at_seg_i_accessor(const int i)
        {
            return segments[i].XNext();
        }
        const XNext_t &XNext_at_seg_i_accessor(const int i) const
        {
            return segments[i].XNext();
        }

        XNextx_t &XNextx_at_seg_i_accessor(const int i)
        {
            return segments[i].XNextx();
        }
        const XNextx_t &XNextx_at_seg_i_accessor(const int i) const
        {
            return segments[i].XNextx();
        }

        XNextw_t &XNextw_at_seg_i_accessor(const int i)
        {
            return segments[i].XNextw();
        }
        const XNextw_t &XNextw_at_seg_i_accessor(const int i) const
        {
            return segments[i].XNextw();
        }

        L_t &L_at_seg_i_accessor(const int i)
        {
            return segments[i].L();
        }
        const L_t &L_at_seg_i_accessor(const int i) const
        {
            return segments[i].L();
        }

        Lx_t &Lx_at_seg_i_accessor(const int i)
        {
            return segments[i].Lx();
        }
        const Lx_t &Lx_at_seg_i_accessor(const int i) const
        {
            return segments[i].Lx();
        }

        Lw_t &Lw_at_seg_i_accessor(const int i)
        {
            return segments[i].Lw();
        }
        const Lw_t &Lw_at_seg_i_accessor(const int i) const
        {
            return segments[i].Lw();
        }

        Lxx_t &Lxx_at_seg_i_accessor(const int i)
        {
            return segments[i].Lxx();
        }
        const Lxx_t &Lxx_at_seg_i_accessor(const int i) const
        {
            return segments[i].Lxx();
        }

        Lxw_t &Lxw_at_seg_i_accessor(const int i)
        {
            return segments[i].Lxw();
        }
        const Lxw_t &Lxw_at_seg_i_accessor(const int i) const
        {
            return segments[i].Lxw();
        }

        Lww_t &Lww_at_seg_i_accessor(const int i)
        {
            return segments[i].Lww();
        }
        const Lww_t &Lww_at_seg_i_accessor(const int i) const
        {
            return segments[i].Lww();
        }

        H_t &H_at_seg_i_accessor(const int i)
        {
            return segments[i].H();
        }
        const H_t &H_at_seg_i_accessor(const int i) const
        {
            return segments[i].H();
        }

        Hx_t &Hx_at_seg_i_accessor(const int i)
        {
            return segments[i].Hx();
        }
        const Hx_t &Hx_at_seg_i_accessor(const int i) const
        {
            return segments[i].Hx();
        }

        Hw_t &Hw_at_seg_i_accessor(const int i)
        {
            return segments[i].Hw();
        }
        const Hw_t &Hw_at_seg_i_accessor(const int i) const
        {
            return segments[i].Hw();
        }

        G_t &G_at_seg_i_accessor(const int i)
        {
            return segments[i].G();
        }
        const G_t &G_at_seg_i_accessor(const int i) const
        {
            return segments[i].G();
        }

        Gx_t &Gx_at_seg_i_accessor(const int i)
        {
            return segments[i].Gx();
        }
        const Gx_t &Gx_at_seg_i_accessor(const int i) const
        {
            return segments[i].Gx();
        }

        Gw_t &Gw_at_seg_i_accessor(const int i)
        {
            return segments[i].Gw();
        }
        const Gw_t &Gw_at_seg_i_accessor(const int i) const
        {
            return segments[i].Gw();
        }

        SegmentDataVector_t segments;
        JumpData_t jump;

    }; // struct PhaseDataDefaultTpl

    template <typename PhaseSpec,
              template <typename PS> class JumpTpl>
    class PhaseModelDefaultTpl
        : public PhaseModelBase<PhaseModelDefaultTpl<PhaseSpec, JumpTpl>, typename PhaseSpec::BS>
    {
    public:
        using PS = PhaseSpec;

        GALILEO_PHASE_SPEC_MASTER_TYPEDEF(PS);

        using Meta_t = PhaseDefaultTpl<PS, JumpTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = PhaseModelBase<PhaseModelDefaultTpl<PS, JumpTpl>, typename PS::BS>;

        using JumpMeta_t = typename traits<Meta_t>::JumpMeta_t;
        using JumpModel_t = typename traits<Meta_t>::JumpModel_t;
        using JumpData_t = typename traits<Meta_t>::JumpData_t;

        PhaseModelDefaultTpl(const PS &ps, const JumpModel_t &jump)
            : Base(),
              ps_(ps),
              jump_(jump)
        {
        }

        template <typename StateMatrixType, typename ControlParamMatrixType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateMatrixType> &xs,
                  const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
        {
            assert(xs.rows() == get_ps().get_nx());
            assert(ws.rows() == get_ps().get_nw());
            assert(ws.cols() == segments_.size());
            assert(xs.cols() == segments_.size() + (with_jump_ ? 1 : 0));

            int i = 0;

            if (with_jump_)
            {
                jump_.calc(data.jump, col(xs, i));
                ++i;
            }

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
            assert(xs.rows() == get_ps().get_nx());
            assert(ws.rows() == get_ps().get_nw());
            assert(ws.cols() == segments_.size());
            assert(xs.cols() == segments_.size() + (with_jump_ ? 1 : 0));

            int i = 0;

            if (with_jump_)
            {
                jump_.calcDiff(data.jump, col(xs, i));
                ++i;
            }
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
            assert(xs.rows() == get_ps().get_nx());
            assert(ws.rows() == get_ps().get_nw());
            assert(ws.cols() == segments_.size());
            assert(xs.cols() == segments_.size() + (with_jump_ ? 1 : 0));

            int i = 0;

            if (with_jump_)
            {
                // null op, but follow same interface as calc/calcDiff
                ++i;
            }

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

        void addSegment(const SegmentModel_t &segment)
        {
            segments_.push_back(segment);
        }

        const PS &get_ps() const
        {
            return ps_.get();
        }

        const SegmentModelVector_t &get_segments() const
        {
            return segments_;
        }

        const JumpModel_t &get_jump() const
        {
            return jump_;
        }

    protected:
        std::reference_wrapper<const PS> ps_;
        SegmentModelVector_t segments_;
        JumpModel_t jump_;

        bool with_jump_ = false;

    }; // class PhaseModelDefaultTpl

} // namespace galileo

#endif // __galileo_predictive_phases_phase_default_hpp__
