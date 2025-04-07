#ifndef __galileo_predictive_phases_phase_generic_hpp__
#define __galileo_predictive_phases_phase_generic_hpp__

#include "galileo/predictive/phases/fwd.hpp"
#include "galileo/predictive/phases/phase-base.hpp"
#include "galileo/predictive/phases/phase-collection.hpp"
#include "galileo/predictive/phases/phase-basic-visitors.hxx"
#include "galileo/utils/aligned-vector.hpp"

#include <boost/mpl/contains.hpp>

namespace galileo
{

    namespace predictive
    {

        template <
            typename PhaseSpec,
            template <typename PS> class PhaseCollectionTpl>
        struct PhaseTpl;

        template <typename PhaseSpec,
                  template <typename PS> class PhaseCollectionTpl>
        struct traits<PhaseTpl<PhaseSpec, PhaseCollectionTpl>>
        {
            using PS = PhaseSpec;
            using PhaseCollection = PhaseCollectionTpl<PS>;

            using PhaseDataDerived = PhaseDataTpl<PS, PhaseCollectionTpl>;
            using PhaseModelDerived = PhaseModelTpl<PS, PhaseCollectionTpl>;

            using SegmentData_t = typename PS::SegmentData_t;
            using SegmentDataVector_t = typename PS::SegmentDataVector_t;
            using SegmentModel_t = typename PS::SegmentModel_t;
        };

        template <typename PhaseSpec,
                  template <typename PS> class PhaseCollectionTpl>
        struct traits<PhaseDataTpl<PhaseSpec, PhaseCollectionTpl>>
        {
            using PS = PhaseSpec;
            using PhaseDerived = PhaseTpl<PS, PhaseCollectionTpl>;
            using PhaseDataDerived = typename traits<PhaseDerived>::PhaseDataDerived;
            using PhaseModelDerived = typename traits<PhaseDerived>::PhaseModelDerived;
        };

        template <typename PhaseSpec,
                  template <typename PS> class PhaseCollectionTpl>
        struct traits<PhaseModelTpl<PhaseSpec, PhaseCollectionTpl>>
        {
            using PS = PhaseSpec;
            using PhaseDerived = PhaseTpl<PS, PhaseCollectionTpl>;
            using PhaseDataDerived = typename traits<PhaseDerived>::PhaseDataDerived;
            using PhaseModelDerived = typename traits<PhaseDerived>::PhaseModelDerived;
        };

        template <typename PhaseSpec,
                  template <typename PS> class PhaseCollectionTpl>
        struct PhaseDataTpl : public PhaseDataBase<PhaseDataTpl<PhaseSpec, PhaseCollectionTpl>, PhaseSpec>,
                              PhaseCollectionTpl<PhaseSpec>::PhaseDataVariant
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using PhaseDerived = PhaseTpl<PS, PhaseCollectionTpl>;
            using PhaseDataDerived = typename traits<PhaseDerived>::PhaseDataDerived;
            using PhaseModelDerived = typename traits<PhaseDerived>::PhaseModelDerived;

            GALILEO_PHASE_DATA_TYPEDEF(PhaseDerived);

            using PhaseCollection = PhaseCollectionTpl<PS>;
            using PhaseDataVariant = typename PhaseCollection::PhaseDataVariant;

            PhaseDataVariant &toVariant()
            {
                return *static_cast<PhaseDataVariant *>(this);
            }
            const PhaseDataVariant &toVariant() const
            {
                return *static_cast<const PhaseDataVariant *>(this);
            }

            SegmentDataVector_t &segments()
            {
                return galileo::predictive::phase_segment_data_vector(*this);
            }

            PhaseDataTpl()
                : PhaseDataVariant()
            {
            }

            PhaseDataTpl(const PhaseDataVariant &phase_data_variant)
                : PhaseDataVariant(phase_data_variant)
            {
            }

            template <typename PhaseDataDerived>
            PhaseDataTpl(const PhaseDataBase<PhaseDataDerived, PhaseSpec> &phase_data)
                : PhaseCollection::PhaseDataVariant((PhaseDataVariant)phase_data.derived())
            {
                BOOST_MPL_ASSERT((boost::mpl::contains<typename PhaseDataVariant::types, PhaseDataDerived>));
            }

            GENERIC_ACCESSOR(SegmentDataVector_t, segments);

        }; // struct PhaseDataTpl

        template <typename PhaseSpec,
                  template <typename PS> class PhaseCollectionTpl>
        struct PhaseModelTpl : public PhaseModelBase<PhaseModelTpl<PhaseSpec, PhaseCollectionTpl>, PhaseSpec>,
                               PhaseCollectionTpl<PhaseSpec>::PhaseModelVariant
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using PS = PhaseSpec;

            using PhaseDerived = PhaseTpl<PS, PhaseCollectionTpl>;
            using PhaseModelDerived = typename traits<PhaseDerived>::PhaseModelDerived;
            using PhaseDataDerived = typename traits<PhaseDerived>::PhaseDataDerived;

            using PhaseCollection = PhaseCollectionTpl<PS>;
            using PhaseModelVariant = typename PhaseCollection::PhaseModelVariant;

            PhaseModelTpl()
                : PhaseModelVariant()
            {
            }

            PhaseModelTpl(const PhaseModelVariant &phase_model_variant)
                : PhaseModelVariant(phase_model_variant)
            {
            }

            template <typename PhaseModelDerived>
            PhaseModelTpl(const PhaseModelBase<PhaseModelDerived, PhaseSpec> &phase_model)
                : PhaseCollection::PhaseModelVariant((PhaseModelVariant)phase_model.derived())
            {
                BOOST_MPL_ASSERT((boost::mpl::contains<typename PhaseModelVariant::types, PhaseModelDerived>));
            }

            PhaseModelVariant &toVariant()
            {
                return *static_cast<PhaseModelVariant *>(this);
            }

            const PhaseModelVariant &toVariant() const
            {
                return *static_cast<const PhaseModelVariant *>(this);
            }

            template <typename StateMatrixType, typename ControlParamMatrixType>
            void calc(PhaseDataDerived &data,
                      const Eigen::MatrixBase<StateMatrixType> &xs,
                      const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
            {
                galileo::predictive::phase_calc_zeroth_order(*this, data, xs.derived(), ws.derived());
            }

            template <typename StateMatrixType, typename ControlParamMatrixType>
            void calcDiff(PhaseDataDerived &data,
                          const Eigen::MatrixBase<StateMatrixType> &xs,
                          const Eigen::MatrixBase<ControlParamMatrixType> &ws) const
            {
                galileo::predictive::phase_calc_first_order(*this, data, xs.derived(), ws.derived());
            }

            template <typename StateMatrixType, typename ControlParamMatrixType>
            void quasiStatic(PhaseDataDerived &data,
                             const Eigen::MatrixBase<StateMatrixType> &xs,
                             Eigen::MatrixBase<ControlParamMatrixType> &ws,
                             const std::size_t maxiter,
                             const typename PS::NumScalar &tol) const
            {
                galileo::predictive::phase_quasi_static(*this, data, xs.derived(), ws.derived(), maxiter, tol);
            }

            template <typename DataCollector>
            auto createData(DataCollector *const collector)
            {
                return galileo::predictive::phase_create_data(*this, collector);
            }

            const typename PS::SegmentModel_t &segment() const
            {
                return galileo::predictive::phase_segment_model(*this);
            }

            const typename PS::NumScalar &period() const
            {
                return galileo::predictive::phase_period(*this);
            }

        }; // struct PhaseModelTpl

    } // namespace predictive

} // namespace galileo

#endif // __galileo_predictive_phases_phase_generic_hpp__
