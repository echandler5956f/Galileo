#include "galileo/fwd.hpp"
#include <Eigen/Dense>

#include "galileo/common/visitors/unary-visitor.hpp"
#include <boost/mpl/contains.hpp>

#include <boost/variant.hpp>

#include <boost/mpl/assert.hpp>

#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <string>

namespace galileo
{

    template <template <typename> class TestMetaTpl, typename VarScalar_>
    struct TestSpecTpl
    {
        using VarScalar = VarScalar_;
        using TS = TestSpecTpl<TestMetaTpl, VarScalar_>;

        using Meta_t = TestMetaTpl<TS>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        static constexpr int NX = traits<Meta_t>::NX;
        static constexpr int NU = traits<Meta_t>::NU;

        using DimNX_t = DimensionTpl<NX>;
        using DimNU_t = DimensionTpl<NU>;

        TestSpecTpl(int nx, int nu)
            : dim_nx(nx), dim_nu(nu)
        {
        }

        DimNX_t dim_nx;
        DimNU_t dim_nu;

        int get_nx() const { return dim_nx.value(); }
        int get_nu() const { return dim_nu.value(); }
    };

    template <typename Derived, typename TestSpec>
    struct TestBaseTpl;

    template <typename Derived, typename TestSpec>
    struct TestDataBase;

    template <typename Derived, typename TestSpec>
    struct TestModelBase;

    template <typename Derived, typename TestSpec>
    struct traits<TestBaseTpl<Derived, TestSpec>>
    {
        using TS = TestSpec;
        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;
    };

    template <typename Derived, typename TestSpec>
    struct traits<TestDataBase<Derived, TestSpec>>
    {
        using TS = TestSpec;
        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;
    };

    template <typename Derived, typename TestSpec>
    struct traits<TestModelBase<Derived, TestSpec>>
    {
        using TS = TestSpec;
        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;
    };

    template <typename Derived, typename TestSpec>
    struct TestDataBase
        : public internal::CRTP<Derived>
    {
    public:
        using TS = TestSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;

        FORWARD_ACCESSOR(MatOut_t, out);

    protected:
        inline TestDataBase()
        {
        }

        inline TestDataBase(const TestDataBase &clone)
        {
        }

        inline TestDataBase &operator=(const TestDataBase &clone)
        {
            return *this;
        }

    }; // struct TestDataBase

    template <typename Derived, typename TestSpec>
    class TestModelBase
        : public internal::CRTP<Derived>
    {
    public:
        using TS = TestSpec;

        using Meta_t = typename traits<Derived>::Meta_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        template <typename StateVectorType, typename ControlVectorType>
            requires IsEigenVector<StateVectorType> && IsEigenVector<ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calc(data, x, u);
        }

        Data_t createData() const
        {
            return this->derived().createData();
        }

        const TS &get_ts() const
        {
            return this->derived().get_ts_impl();
        }

        const TS &get_ts_impl() const
        {
            return ts_.get();
        }

    protected:
        inline TestModelBase(const TS &ts)
            : ts_(ts)
        {
        }

        inline TestModelBase(const TestModelBase &clone)
            : ts_(clone.ts_)
        {
        }

        inline TestModelBase &operator=(const TestModelBase &clone)
        {
            ts_ = clone.ts_;
            return *this;
        }

        std::reference_wrapper<const TS> ts_;

    }; // class TestModelBase

    template <typename TestSpec>
    struct TestDerivedTpl;

    template <typename TestSpec>
    struct TestDerivedDataTpl;

    template <typename TestSpec>
    struct TestDerivedModelTpl;

    template <typename TestSpec>
    struct traits<TestDerivedTpl<TestSpec>>
    {
        using TS = TestSpec;
        using SpecOfBaseClass = TS;
        using Meta_t = TestDerivedTpl<TS>;
        using Data_t = TestDerivedDataTpl<TS>;
        using Model_t = TestDerivedModelTpl<TS>;

        using VarScalar = typename TS::VarScalar;
        static constexpr int NX = 5;
        static constexpr int NU = 3;

        using MatOut_t = Eigen::GMatrix<VarScalar, NX, NU>;
    };

    template <typename TestSpec>
    struct traits<TestDerivedDataTpl<TestSpec>>
    {
        using TS = TestSpec;
        using SpecOfBaseClass = TS;
        using Meta_t = TestDerivedTpl<TS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;
    };

    template <typename TestSpec>
    struct traits<TestDerivedModelTpl<TestSpec>>
    {
        using TS = TestSpec;
        using SpecOfBaseClass = TS;
        using Meta_t = TestDerivedTpl<TS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;
    };

    template <typename TestSpec>
    struct TestDerivedDataTpl
        : public TestDataBase<TestDerivedDataTpl<TestSpec>, TestSpec>
    {
    public:
        using TS = TestSpec;

        using Meta_t = TestDerivedTpl<TS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Base = TestDataBase<TestDerivedDataTpl<TS>, TS>;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;

        DEFAULT_ACCESSOR(MatOut_t, out);

        TestDerivedDataTpl(const Model_t &model)
            : Base(),
              out(model.get_ts().get_nx(), model.get_ts().get_nu())
        {
            out.setZero();
        }

        MatOut_t out;
    };

    template <typename TestSpec>
    struct TestDerivedModelTpl
        : public TestModelBase<TestDerivedModelTpl<TestSpec>, TestSpec>
    {
    public:
        using TS = TestSpec;

        using Meta_t = TestDerivedTpl<TS>;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Base = TestModelBase<TestDerivedModelTpl<TS>, TS>;

        using VarScalar = typename TS::VarScalar;
        using DimNX_t = typename TS::DimNX_t;
        using DimNU_t = typename TS::DimNU_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;

        TestDerivedModelTpl(const TS &ts, VarScalar some_constant)
            : Base(ts),
              some_constant_(some_constant)
        {
        }

        // dummy operation for testing
        template <typename StateVectorType, typename ControlVectorType>
            requires IsEigenVector<StateVectorType> && IsEigenVector<ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            static_assert(StateVectorType::SizeAtCompileTime == DimNX_t::Value, "State vector size mismatch");
            static_assert(ControlVectorType::SizeAtCompileTime == DimNU_t::Value, "Control vector size mismatch");
            static_assert(std::is_same<VarScalar, typename StateVectorType::Scalar>::value, "State vector scalar mismatch");
            static_assert(std::is_same<VarScalar, typename ControlVectorType::Scalar>::value, "Control vector scalar mismatch");

            assert(get_ts().get_nx() == x.size());
            assert(get_ts().get_nu() == u.size());

            MatOut_t M = x * u.transpose() + MatOut_t::Identity() * some_constant_;
            auto A1 = M.array().sin();
            auto A2 = (x.replicate(1, get_ts().get_nu()).array() + u.transpose().replicate(get_ts().get_nx(), 1).array()).cos();
            data.out = (A1 + A2).matrix();
        }

        Data_t createData() const
        {
            return Data_t(*this);
        }

        using Base::get_ts;

    protected:
        VarScalar some_constant_;
    };

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct TestDataTpl;

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct TestModelTpl;

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct TestTpl;

    namespace fusion
    {
        struct TestFamily
        {
        };

        template <>
        struct UnaryVisitorFamilyTraits<TestFamily>
        {
            template <typename TS, template <typename> class CollectionTpl>
            using ModelTpl = TestModelTpl<TS, CollectionTpl>;

            template <typename TS, template <typename> class CollectionTpl>
            using DataTpl = TestDataTpl<TS, CollectionTpl>;

            template <typename ModelType, typename TS>
            using ModelBase = TestModelBase<ModelType, TS>;

            template <typename DataType, typename TS>
            using DataBase = TestDataBase<DataType, TS>;
        };

        template <typename TestVisitorDerived, typename ReturnType = void>
        using TestUnaryVisitorBase = UnaryVisitorBase<TestFamily, TestVisitorDerived, ReturnType>;
    } // namespace fusion

    // --------------------------------------------------------------------------------------------------------

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    inline void test_calc(
        const TestModelTpl<TestSpec, TestCollectionTpl> &test_model,
        TestDataTpl<TestSpec, TestCollectionTpl> &test_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Eigen::MatrixBase<ControlVectorType> &u);

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    struct TestCalcVisitor
        : fusion::TestUnaryVisitorBase<TestCalcVisitor<TestSpec, TestCollectionTpl, StateVectorType, ControlVectorType>>
    {
        using ArgsType = boost::fusion::vector<const StateVectorType &, const ControlVectorType &>;

        template <typename TestModelType>
        static void algo(
            const TestModelBase<TestModelType, TestSpec> &test_model,
            TestDataBase<typename traits<TestModelType>::Data_t, TestSpec> &test_data,
            const Eigen::MatrixBase<StateVectorType> &x,
            const Eigen::MatrixBase<ControlVectorType> &u)
        {
            test_model.calc(test_data.derived(), x.derived(), u.derived());
        }
    };

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl,
              typename StateVectorType,
              typename ControlVectorType>
    inline void test_calc(
        const TestModelTpl<TestSpec, TestCollectionTpl> &test_model,
        TestDataTpl<TestSpec, TestCollectionTpl> &test_data,
        const Eigen::MatrixBase<StateVectorType> &x,
        const Eigen::MatrixBase<ControlVectorType> &u)
    {
        typedef TestCalcVisitor<TestSpec, TestCollectionTpl, StateVectorType, ControlVectorType> Algo;

        Algo::run(test_model, test_data, typename Algo::ArgsType(x.derived(), u.derived()));
    }

    // --------------------------------------------------------------------------------------------------------

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    inline TestDataTpl<TestSpec, TestCollectionTpl> test_create_data(
        const TestModelTpl<TestSpec, TestCollectionTpl> &test_model);

    template <typename TestSpec,
              template <typename> class TestCollectionTpl>
    struct TestCreateDataVisitor
        : fusion::TestUnaryVisitorBase<TestCreateDataVisitor<TestSpec, TestCollectionTpl>,
                                       TestDataTpl<TestSpec, TestCollectionTpl>>
    {
        using TestCollection_t = TestCollectionTpl<TestSpec>;
        using TestModelVariant_t = TestCollection_t::TestModelVariant_t;
        using TestDataVariant_t = TestDataTpl<TestSpec, TestCollectionTpl>;

        template <typename TestModelType>
        static TestDataVariant_t algo(
            const TestModelBase<TestModelType, TestSpec> &test_model)
        {
            return TestDataVariant_t(test_model.createData());
        }
    };

    template <typename TestSpec,
              template <typename> class TestCollectionTpl>
    inline TestDataTpl<TestSpec, TestCollectionTpl> test_create_data(
        const TestModelTpl<TestSpec, TestCollectionTpl> &test_model)
    {
        typedef TestCreateDataVisitor<TestSpec, TestCollectionTpl> Algo;
        return Algo::run(test_model);
    }

    // --------------------------------------------------------------------------------------------------------

    template <typename TestSpec, template <typename> class TestCollectionTpl>
    inline TestSpec test_get_ts(const TestModelTpl<TestSpec, TestCollectionTpl> &test_model);

    template <typename TestSpec, template <typename> class TestCollectionTpl>
    struct TestGetTestSpecVisitor
        : boost::static_visitor<const TestSpec &>
    {

        using ReturnType = const TestSpec &;

        template <typename TestModelType>
        ReturnType operator()(const TestModelBase<TestModelType, TestSpec> &test_model) const
        {
            return test_model.get_ts();
        }

        static ReturnType run(const TestModelTpl<TestSpec, TestCollectionTpl> &test_model)
        {
            return boost::apply_visitor(TestGetTestSpecVisitor<TestSpec, TestCollectionTpl>(), test_model);
        }
    };

    template <typename TestSpec, template <typename> class TestCollectionTpl>
    inline const TestSpec &test_get_ts(const TestModelTpl<TestSpec, TestCollectionTpl> &test_model)
    {
        return TestGetTestSpecVisitor<TestSpec, TestCollectionTpl>::run(test_model);
    }

    // --------------------------------------------------------------------------------------------------------

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    inline typename TestDataTpl<TestSpec, TestCollectionTpl>::MatOut_t test_out(
        const TestDataTpl<TestSpec, TestCollectionTpl> &test_data);

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct TestOutVisitor : boost::static_visitor<typename TestDataTpl<TestSpec, TestCollectionTpl>::MatOut_t>
    {

        using ReturnType = typename TestDataTpl<TestSpec, TestCollectionTpl>::MatOut_t;

        template <typename TestDataType>
        ReturnType operator()(const TestDataBase<TestDataType, TestSpec> &test_data) const
        {
            return test_data.out();
        }

        static ReturnType run(const TestDataTpl<TestSpec, TestCollectionTpl> &test_data)
        {
            return boost::apply_visitor(TestOutVisitor<TestSpec, TestCollectionTpl>(), test_data);
        }
    };

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    inline typename TestDataTpl<TestSpec, TestCollectionTpl>::MatOut_t test_out(
        const TestDataTpl<TestSpec, TestCollectionTpl> &test_data)
    {
        return TestOutVisitor<TestSpec, TestCollectionTpl>::run(test_data);
    }

    // --------------------------------------------------------------------------------------------------------

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct traits<TestTpl<TestSpec, TestCollectionTpl>>
    {
        using TS = TestSpec;
        using SpecOfBaseClass = TS;

        using Meta_t = TestTpl<TS, TestCollectionTpl>;
        using Collection_t = TestCollectionTpl<TS>;
        using Model_t = TestModelTpl<TS, TestCollectionTpl>;
        using Data_t = TestDataTpl<TS, TestCollectionTpl>;

        using VarScalar = typename TS::VarScalar;

        static constexpr int NX = Eigen::Dynamic;
        static constexpr int NU = Eigen::Dynamic;

        using MatOut_t = Eigen::GMatrix<VarScalar, NX, NU>;
    };

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct traits<TestDataTpl<TestSpec, TestCollectionTpl>>
    {
        using TS = TestSpec;
        using SpecOfBaseClass = TS;
        using Meta_t = TestTpl<TS, TestCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;
    };

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct traits<TestModelTpl<TestSpec, TestCollectionTpl>>
    {
        using TS = TestSpec;
        using SpecOfBaseClass = TS;
        using Meta_t = TestTpl<TS, TestCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Model_t = typename traits<Meta_t>::Model_t;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;
    };

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct TestDataTpl
        : public TestDataBase<TestDataTpl<TestSpec, TestCollectionTpl>, TestSpec>,
          TestCollectionTpl<TestSpec>::TestDataVariant_t
    {
    public:
        using TS = TestSpec;
        using Meta_t = TestTpl<TS, TestCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = TestDataBase<TestDataTpl<TS, TestCollectionTpl>, TS>;

        using MatOut_t = typename traits<Meta_t>::MatOut_t;

        using DataVariant_t = typename Collection_t::TestDataVariant_t;

        DataVariant_t &toVariant()
        {
            return *static_cast<DataVariant_t *>(this);
        }
        const DataVariant_t &toVariant() const
        {
            return *static_cast<const DataVariant_t *>(this);
        }

        MatOut_t out() const
        {
            return galileo::test_out(*this);
        }

        TestDataTpl()
            : DataVariant_t()
        {
        }

        TestDataTpl(const DataVariant_t &data_variant)
            : DataVariant_t(data_variant)
        {
        }

        template <typename DataDerived>
        TestDataTpl(const TestDataBase<DataDerived, TestSpec> &data)
            : Collection_t::TestDataVariant_t((DataVariant_t)data.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename DataVariant_t::types, DataDerived>));
        }

        GENERIC_ACCESSOR(MatOut_t, out);
    };

    template <typename TestSpec,
              template <typename TS> class TestCollectionTpl>
    struct TestModelTpl
        : public TestModelBase<TestModelTpl<TestSpec, TestCollectionTpl>, TestSpec>,
          TestCollectionTpl<TestSpec>::TestModelVariant_t
    {
    public:
        using TS = TestSpec;
        using Meta_t = TestTpl<TS, TestCollectionTpl>;
        using Collection_t = typename traits<Meta_t>::Collection_t;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;
        using Base = TestModelBase<TestModelTpl<TS, TestCollectionTpl>, TS>;

        using ModelVariant_t = typename Collection_t::TestModelVariant_t;

        ModelVariant_t &toVariant()
        {
            return *static_cast<ModelVariant_t *>(this);
        }

        const ModelVariant_t &toVariant() const
        {
            return *static_cast<const ModelVariant_t *>(this);
        }

        TestModelTpl()
            : ModelVariant_t()
        {
        }

        TestModelTpl(const ModelVariant_t &model_variant)
            : ModelVariant_t(model_variant)
        {
        }

        template <typename ModelDerived>
        TestModelTpl(const TestModelBase<ModelDerived, TestSpec> &model)
            : Base(model.get_ts()),
              ModelVariant_t((ModelVariant_t)model.derived())
        {
            BOOST_MPL_ASSERT((boost::mpl::contains<typename ModelVariant_t::types, ModelDerived>));
        }

        template <typename StateVectorType, typename ControlVectorType>
            requires IsEigenVector<StateVectorType> && IsEigenVector<ControlVectorType>
        void calc(Data_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            galileo::test_calc(*this, data, x.derived(), u.derived());
        }

        Data_t createData() const
        {
            return galileo::test_create_data(*this);
        }

        using Base::get_ts;

        const TS &get_ts_impl() const
        {
            return galileo::test_get_ts(*this);
        }
    };

} // namespace galileo

template <typename TestSpec>
using DerivedMetaTpl = galileo::TestDerivedTpl<TestSpec>;
template <typename TestSpec>
using DerivedModelTpl = galileo::TestDerivedModelTpl<TestSpec>;
template <typename TestSpec>
using DerivedDataTpl = galileo::TestDerivedDataTpl<TestSpec>;

template <typename TestSpec>
struct TestCollectionDefaultTpl
{
    using TestModelVariant_t = boost::variant<DerivedModelTpl<TestSpec>>;
    using TestDataVariant_t = boost::variant<DerivedDataTpl<TestSpec>>;
}; // struct TestCollectionDefaultTpl

using VarScalar = double;
using TestSpec_t = galileo::TestSpecTpl<DerivedMetaTpl, VarScalar>;
using TS = TestSpec_t;

using Model_t = typename TS::Model_t;
using Data_t = typename TS::Data_t;

using VectorNx_t = Eigen::GMatrix<VarScalar, TS::NX, 1>;
using VectorNu_t = Eigen::GMatrix<VarScalar, TS::NU, 1>;

using ModelTpl = galileo::TestModelTpl<TS, TestCollectionDefaultTpl>;
using DataTpl = galileo::TestDataTpl<TS, TestCollectionDefaultTpl>;

int main()
{
    TS ts(5, 3);

    std::vector<ModelTpl> models;
    std::vector<DataTpl> datas;

    std::uniform_real_distribution<VarScalar> unif(0.0, 1.0);
    std::default_random_engine re;

    for (int i = 0; i < 10; ++i)
    {
        VarScalar some_constant = unif(re);
        models.push_back(Model_t(ts, some_constant));
        datas.push_back(models.back().createData());
    }

    VectorNx_t x;
    VectorNu_t u;

    for (int i = 0; i < 10; ++i)
    {
        x = VectorNx_t::Random(ts.get_nx());
        u = VectorNu_t::Random(ts.get_nu());

        models[i].calc(datas[i], x, u);
    }

    for (int i = 0; i < 10; ++i)
    {
        std::cout << datas[i].out() << std::endl;
    }

    return 0;
}
