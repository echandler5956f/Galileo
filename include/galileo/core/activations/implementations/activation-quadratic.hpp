#ifndef __galileo_core_activations_activation_quadratic_hpp__
#define __galileo_core_activations_activation_quadratic_hpp__

#include "galileo/core/activations/activation-base.hpp"

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ResidualTpl>
    struct ActivationQuadraticTpl;

    template <typename PhaseSpec,
              template <typename PS> class ResidualTpl>
    struct traits<ActivationQuadraticTpl<PhaseSpec, ResidualTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ActivationQuadraticTpl<PS, ResidualTpl>;
        using Model_t = ActivationModelQuadraticTpl<PS, ResidualTpl>;
        using Data_t = ActivationDataQuadraticTpl<PS, ResidualTpl>;

        using ResidualMeta_t = typename traits<ResidualTpl<PS>>::Meta_t;
        using ResidualModel_t = typename traits<ResidualMeta_t>::Model_t;
        using ResidualData_t = typename traits<ResidualMeta_t>::Data_t;

        using DimNR_t = typename traits<ResidualMeta_t>::DimNR_t;

        using A_t = typename PS::VarScalar;
        using Ar_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, 1, PS::Options>;
        using Arr_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, DimNR_t::Value, PS::Options>;
        using Arr_diag_t = Eigen::DiagonalMatrix<typename PS::VarScalar, DimNR_t::Value>;
    };

    template <typename PhaseSpec,
              template <typename PS> class ResidualTpl>
    struct traits<ActivationDataQuadraticTpl<PhaseSpec, ResidualTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ActivationQuadraticTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ResidualTpl>
    struct traits<ActivationModelQuadraticTpl<PhaseSpec, ResidualTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ActivationQuadraticTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ResidualTpl>
    struct ActivationDataQuadraticTpl
        : public ActivationDataBase<ActivationDataQuadraticTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ActivationQuadraticTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActivationDataBase<ActivationDataQuadraticTpl<PS, ResidualTpl>, PS>;

        GALILEO_ACTIVATION_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(A_t, A);
        DEFAULT_ACCESSOR(Ar_t, Ar);
        DEFAULT_ACCESSOR(Arr_t, Arr);

        ActivationDataQuadraticTpl(const Model_t &model)
            : A(0.), Ar(model.get_nr()), Arr(Arr_diag_t(model.get_nr()))
        {
            Ar.setZero();
            Arr.setIdentity();
        }

        A_t A;
        Ar_t Ar;
        Arr_t Arr;

    }; // class ActivationDataQuadraticTpl

    template <typename PhaseSpec,
              template <typename PS> class ResidualTpl>
    class ActivationModelQuadraticTpl
        : public ActivationModelBase<ActivationModelQuadraticTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ActivationQuadraticTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActivationModelBase<ActivationModelQuadraticTpl<PS, ResidualTpl>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;

        explicit ActivationModelQuadraticTpl(const PS &ps, const DimNR_t &nr_dim)
            : Base(ps, nr_dim)
        {
        }

        template <typename ResidualVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            std::cout << "activation quadratic calc" << std::endl;
            data.A = typename PS::VarScalar(0.5) * r.dot(r);
            std::cout << "data.A: " << data.A << std::endl;
        }

        template <typename ResidualVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            data.Ar = r;
            // The Hessian has constant values which were set in createData.
        }

        Data_t createData() const
        {
            return Data_t(*this);
        }

        using Base::get_ps;

        using Base::get_nr;
        using Base::get_nr_dim;

    }; // class ActivationModelQuadraticTpl

} // namespace galileo

#endif // __galileo_core_activations_activation_quadratic_hpp__
