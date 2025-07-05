#ifndef __galileo_core_constraints_constraint_manager_hpp__
#define __galileo_core_constraints_constraint_manager_hpp__

#include <iostream>
#include <map>
#include <set>
#include <string>

#include "galileo/core/constraints/fwd.hpp"
#include "galileo/core/constraints/constraint-generic.hpp"

// Despite its name, this file only pertains to constraint residuals,
// not just any general constraints derived from ConstraintModelBase/ConstraintDataBase

namespace galileo
{

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintManagerTpl;

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct ConstraintItemTpl
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        ConstraintItemTpl() {}
        ConstraintItemTpl(const std::string &name_in, const Model_t &model_in, bool active_in = true)
            : name(name_in), model(model_in), active(active_in) {}

        std::string name;
        Model_t model;
        bool active;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintManagerTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = ConstraintCollectionTpl<PS>;
        using ModelManager_t = ConstraintModelManagerTpl<PS, ConstraintCollectionTpl>;
        using DataManager_t = ConstraintDataManagerTpl<PS, ConstraintCollectionTpl>;

        using Meta_t = ConstraintTpl<PS, ConstraintCollectionTpl>;
        using Model_t = typename traits<Meta_t>::Model_t;
        using Data_t = typename traits<Meta_t>::Data_t;

        using Item_t = ConstraintItemTpl<PS, ConstraintCollectionTpl>;

        using ModelContainer_t = std::map<std::string, Item_t>;
        using DataContainer_t = std::map<std::string, Data_t>;

        static constexpr int NH = Eigen::Dynamic;
        static constexpr int NG = Eigen::Dynamic;

        using H_t = Eigen::GMatrix<typename PS::VarScalar, NH, 1, PS::Options>;
        using Hx_t = Eigen::GMatrix<typename PS::VarScalar, NH, PS::NDX, PS::Options>;
        using Hu_t = Eigen::GMatrix<typename PS::VarScalar, NH, PS::NU, PS::Options>;
        using G_t = Eigen::GMatrix<typename PS::VarScalar, NG, 1, PS::Options>;
        using Gx_t = Eigen::GMatrix<typename PS::VarScalar, NG, PS::NDX, PS::Options>;
        using Gu_t = Eigen::GMatrix<typename PS::VarScalar, NG, PS::NU, PS::Options>;

        using BoundVector_t = Eigen::GMatrix<typename PS::NumScalar, NG, 1, PS::Options>;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintDataManagerTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    struct traits<ConstraintModelManagerTpl<PhaseSpec, ConstraintCollectionTpl>>
    {
        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;
    };

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    class ConstraintDataManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        GALILEO_CONSTRAINT_DATA_TYPEDEF(MetaManager_t);

        DataContainer_t constraints;
        Eigen::Map<H_t> H;
        Eigen::Map<Hx_t> Hx;
        Eigen::Map<Hu_t> Hu;
        Eigen::Map<G_t> G;
        Eigen::Map<Gx_t> Gx;
        Eigen::Map<Gu_t> Gu;

    }; // class ConstraintDataManagerTpl

    template <typename PhaseSpec,
              template <typename PS> class ConstraintCollectionTpl>
    class ConstraintModelManagerTpl
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using PS = PhaseSpec;

        using MetaManager_t = ConstraintManagerTpl<PS, ConstraintCollectionTpl>;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        using BoundVector_t = typename traits<MetaManager_t>::BoundVector_t;

        ConstraintModelManagerTpl() : nh_(0), ng_(0) {}

        void addConstraint(const std::string &name, const Model_t &model, const bool active = true)
        {
            std::pair<typename ModelContainer_t::iterator, bool> ret =
                constraints_.insert(std::make_pair(
                    name, Item_t(name, model, active)));
            if (ret.second == false)
            {
                std::cout << "Warning: we couldn't add the " << name
                          << " constraint item, it already existed." << std::endl;
            }
            else if (active)
            {
                ng_ += model.ng();
                nh_ += model.nh();
                active_set_.insert(name);
                lb_.resize(ng_);
                ub_.resize(ng_);
            }
            else if (!active)
            {
                inactive_set_.insert(name);
            }
        }

        void removeConstraint(const std::string &name)
        {
            typename ModelContainer_t::iterator it = constraints_.find(name);
            if (it != constraints_.end())
            {
                ng_ -= it->second.model.ng();
                nh_ -= it->second.model.nh();
                constraints_.erase(it);
                inactive_set_.erase(name);
                lb_.resize(ng_);
                ub_.resize(ng_);
            }
            else
            {
                std::cout << "Warning: we couldn't remove the " << name
                          << " constraint item, it doesn't exist." << std::endl;
            }
        }

        void changeConstraintStatus(const std::string &name, bool active)
        {
            typename ModelContainer_t::iterator it = constraints_.find(name);
            if (it != constraints_.end())
            {
                if (active && !it->second.active)
                {
                    ng_ += it->second.model.ng();
                    nh_ += it->second.model.nh();
                    active_set_.insert(name);
                    inactive_set_.erase(name);
                    it->second.active = active;
                    lb_.resize(ng_);
                    ub_.resize(ng_);
                }
                else if (!active && it->second.active)
                {
                    ng_ -= it->second.model.ng();
                    nh_ -= it->second.model.nh();
                    active_set_.erase(name);
                    inactive_set_.insert(name);
                    it->second.active = active;
                    lb_.resize(ng_);
                    ub_.resize(ng_);
                }
            }
            else
            {
                std::cout << "Warning: we couldn't change the status of the " << name
                          << " constraint item, it doesn't exist." << std::endl;
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u)
        {
            std::size_t ng_i = 0;
            std::size_t nh_i = 0;

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = constraints_.begin(), end_m = constraints_.end(),
                it_d = data.constraints.begin(), end_d = data.constraints.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x.derived(), u.derived());
                    const std::size_t ng = m_i.model.ng();
                    const std::size_t nh = m_i.model.nh();
                    data.G.segment(ng_i, ng) = d_i.G();
                    data.H.segment(nh_i, nh) = d_i.H();
                    lb_.segment(ng_i, ng) = m_i.model.lb();
                    ub_.segment(ng_i, ng) = m_i.model.ub();
                    ng_i += ng;
                    nh_i += nh;
                }
            }
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x)
        {
            std::size_t ng_i = 0;
            std::size_t nh_i = 0;

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = constraints_.begin(), end_m = constraints_.end(),
                it_d = data.constraints.begin(), end_d = data.constraints.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calc(d_i, x.derived());
                    const std::size_t ng = m_i.model.ng();
                    const std::size_t nh = m_i.model.nh();
                    data.G.segment(ng_i, ng) = d_i.G();
                    data.H.segment(nh_i, nh) = d_i.H();
                    lb_.segment(ng_i, ng) = m_i.model.lb();
                    ub_.segment(ng_i, ng) = m_i.model.ub();
                    ng_i += ng;
                    nh_i += nh;
                }
            }
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u)
        {
            std::size_t ng_i = 0;
            std::size_t nh_i = 0;

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = constraints_.begin(), end_m = constraints_.end(),
                it_d = data.constraints.begin(), end_d = data.constraints.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x.derived(), u.derived());
                    const std::size_t ng = m_i.model.ng();
                    const std::size_t nh = m_i.model.nh();
                    data.Gx.block(ng_i, 0, ng, PS::NDX) = d_i.Gx();
                    data.Gu.block(ng_i, 0, ng, PS::NU) = d_i.Gu();
                    data.Hx.block(nh_i, 0, nh, PS::NDX) = d_i.Hx();
                    data.Hu.block(nh_i, 0, nh, PS::NU) = d_i.Hu();
                    ng_i += ng;
                    nh_i += nh;
                }
            }
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x)
        {
            std::size_t ng_i = 0;
            std::size_t nh_i = 0;

            typename ModelContainer_t::iterator it_m, end_m;
            typename DataContainer_t::iterator it_d, end_d;
            for (it_m = constraints_.begin(), end_m = constraints_.end(),
                it_d = data.constraints.begin(), end_d = data.constraints.end();
                 it_m != end_m || it_d != end_d; ++it_m, ++it_d)
            {
                const Item_t &m_i = it_m->second;
                if (m_i.active)
                {
                    Data_t &d_i = it_d->second;

                    m_i.model.calcDiff(d_i, x.derived());
                    const std::size_t ng = m_i.model.ng();
                    const std::size_t nh = m_i.model.nh();
                    data.Gx.block(ng_i, 0, ng, PS::NDX) = d_i.Gx();
                    data.Hx.block(nh_i, 0, nh, PS::NDX) = d_i.Hx();
                    ng_i += ng;
                    nh_i += nh;
                }
            }
        }

        const std::set<std::string> &getActiveSet() const
        {
            return active_set_;
        }

        const std::set<std::string> &getInactiveSet() const
        {
            return inactive_set_;
        }

        const BoundVector_t &lb() const
        {
            return lb_;
        }

        const BoundVector_t &ub() const
        {
            return ub_;
        }

        bool getConstraintStatus(const std::string &name) const
        {
            typename ModelContainer_t::const_iterator it =
                constraints_.find(name);
            if (it != constraints_.end())
            {
                return it->second.active;
            }
            else
            {
                std::cout << "Warning: we couldn't get the status of the " << name
                          << " constraint item, it doesn't exist." << std::endl;
                return false;
            }
        }

    protected:
        ModelContainer_t constraints_;
        BoundVector_t lb_;
        BoundVector_t ub_;

        std::size_t nh_;
        std::size_t ng_;

        std::set<std::string> active_set_;
        std::set<std::string> inactive_set_;

    }; // class ConstraintModelManagerTpl

} // namespace galileo

#endif // __galileo_core_constraints_constraint_manager_hpp__
