#ifndef __galileo_core_activations_quadratic_hpp__
#define __galileo_core_activations_quadratic_hpp__

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

        static constexpr int NR = traits<ResidualMeta_t>::NR;

        using A_t = PS::VarScalar;
        using Ar_t = Eigen::Matrix<typename PS::VarScalar, NR, 1, PS::Options>;
        using Arr_t = Eigen::Matrix<typename PS::VarScalar, NR, NR, PS::Options>;
        using Arr_diag_t = Eigen::DiagonalMatrix<typename PS::VarScalar, NR>;
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
    struct ActivationDataQuadraticTpl : public ActivationDataBase<ActivationDataQuadraticTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ActivationQuadraticTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        GALILEO_ACTIVATION_DATA_TYPEDEF(Meta_t);

        DEFAULT_ACCESSOR(A_t, A);
        DEFAULT_ACCESSOR(Ar_t, Ar);
        DEFAULT_ACCESSOR(Arr_t, Arr);

        ActivationDataQuadraticTpl() : A(A_t(0.)), Ar(Ar_t::Zero()), Arr(Arr_diag_t())
        {
            Arr.setZero();
        }

        A_t A;
        Ar_t Ar;
        Arr_t Arr;

    }; // class ActivationDataQuadraticTpl

    template <typename PhaseSpec,
              template <typename PS> class ResidualTpl>
    class ActivationModelQuadraticTpl : public ActivationModelBase<ActivationModelQuadraticTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using Meta_t = ActivationQuadraticTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        explicit ActivationModelQuadraticTpl(const int nr) : nr_(nr) {}

        template <typename ResidualVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            data.A() = typename PS::VarScalar(0.5) * r.dot(r);
        }

        template <typename ResidualVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            data.Ar() = r;
            // The Hessian has constant values which were set in createData.
        }

        Data_t createData() const
        {
            Data_t data = Data_t();
            data.Arr.diagonal().setOnes();
            return data;
        }

        int nr_impl() const
        {
            return nr_;
        }

    protected:
        int nr_;

    }; // class ActivationModelQuadraticTpl

} // namespace galileo

#endif // __galileo_core_activations_quadratic_hpp__