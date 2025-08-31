#ifndef __galileo_core_activations_activation_weighted_quadratic_hpp__
#define __galileo_core_activations_activation_weighted_quadratic_hpp__

#include "galileo/core/activations/activation-base.hpp"

namespace galileo
{

    template <typename PhaseSpec, template <typename> class ResidualTpl>
    struct ActivationWeightedQuadraticTpl;

    template <typename PhaseSpec, template <typename> class ResidualTpl>
    struct traits<ActivationWeightedQuadraticTpl<PhaseSpec, ResidualTpl>>
    {
        using PS = PhaseSpec;

        using Meta_t = ActivationWeightedQuadraticTpl<PS, ResidualTpl>;
        using Model_t = ActivationModelWeightedQuadraticTpl<PS, ResidualTpl>;
        using Data_t = ActivationDataWeightedQuadraticTpl<PS, ResidualTpl>;

        using ResidualMeta_t = typename traits<ResidualTpl<PS>>::Meta_t;
        using ResidualModel_t = typename traits<ResidualMeta_t>::Model_t;
        using ResidualData_t = typename traits<ResidualMeta_t>::Data_t;

        using DimNR_t = typename traits<ResidualMeta_t>::DimNR_t;

        using A_t = typename PS::VarScalar;
        using Ar_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, 1, PS::Options>>;
        using Arr_t = ArenaMatrixTpl<Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, DimNR_t::Value, PS::Options>>;

        using WeightVector_t = Eigen::GMatrix<typename PS::VarScalar, DimNR_t::Value, 1, PS::Options>;
    };

    template <typename PhaseSpec, template <typename> class ResidualTpl>
    struct traits<ActivationDataWeightedQuadraticTpl<PhaseSpec, ResidualTpl>>
    {
        using Meta_t = ActivationWeightedQuadraticTpl<PhaseSpec, ResidualTpl>;
    };

    template <typename PhaseSpec, template <typename> class ResidualTpl>
    struct traits<ActivationModelWeightedQuadraticTpl<PhaseSpec, ResidualTpl>>
    {
        using Meta_t = ActivationWeightedQuadraticTpl<PhaseSpec, ResidualTpl>;
    };

    template <typename PhaseSpec, template <typename> class ResidualTpl>
    struct ActivationDataWeightedQuadraticTpl
        : public ActivationDataBase<ActivationDataWeightedQuadraticTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ActivationWeightedQuadraticTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActivationDataBase<ActivationDataWeightedQuadraticTpl<PS, ResidualTpl>, PS>;

        GALILEO_ACTIVATION_DATA_TYPEDEF(Meta_t);

        using WeightVector_t = typename traits<Meta_t>::WeightVector_t;

        DEFAULT_ACCESSOR(A_t, A);
        DEFAULT_ACCESSOR(Ar_t, Ar);
        DEFAULT_ACCESSOR(Arr_t, Arr);

        ActivationDataWeightedQuadraticTpl(const Model_t &model, MemoryArena &arena)
            : A(0.),
              Ar(arena, model.get_nr(), 1),
              Wr(arena, model.get_nr(), 1),
              Arr(arena, model.get_nr(), model.get_nr())
        {
            Ar.setZero();
            Wr.setZero();
            Arr.diagonal() = model.get_weights();
        }

        A_t A;
        Ar_t Ar;
        ArenaMatrixTpl<WeightVector_t> Wr;
        Arr_t Arr;

    }; // class ActivationDataWeightedQuadraticTpl

    template <typename PhaseSpec, template <typename> class ResidualTpl>
    class ActivationModelWeightedQuadraticTpl
        : public ActivationModelBase<ActivationModelWeightedQuadraticTpl<PhaseSpec, ResidualTpl>, PhaseSpec>
    {
    public:
        using PS = PhaseSpec;

        using Meta_t = ActivationWeightedQuadraticTpl<PS, ResidualTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = ActivationModelBase<ActivationModelWeightedQuadraticTpl<PS, ResidualTpl>, PS>;

        using DimNR_t = typename traits<Meta_t>::DimNR_t;
        using WeightVector_t = typename traits<Meta_t>::WeightVector_t;

        explicit ActivationModelWeightedQuadraticTpl(const PS &ps, const DimNR_t &nr_dim, const WeightVector_t &weights)
            : Base(nr_dim), weights_(weights)
        {
        }

        template <typename ResidualVectorType>
        void calc(Data_t &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            data.Wr = weights_.cwiseProduct(r);
            data.A = typename PS::VarScalar(0.5) * r.dot(data.Wr);
        }

        template <typename ResidualVectorType>
        void calcDiff(Data_t &data, const Eigen::MatrixBase<ResidualVectorType> &r) const
        {
            data.Ar = data.Wr;
            data.Arr.diagonal() = weights_;
        }

        Data_t createData(MemoryArena &arena) const { return Data_t(*this, arena); }

        const WeightVector_t &get_weights() const { return weights_; }
        void setWeights(const WeightVector_t &weights) { weights_ = weights; }

        using Base::get_nr;
        using Base::get_nr_dim;

    protected:
        WeightVector_t weights_;

    }; // class ActivationModelWeightedQuadraticTpl

} // namespace galileo

#endif // __galileo_core_activations_activation_weighted_quadratic_hpp__
