#ifndef __galileo_common_container_manager_base_hpp__
#define __galileo_common_container_manager_base_hpp__

#include "galileo/fwd.hpp"

#include <iostream>
#include <map>
#include <set>
#include <string>

namespace galileo
{

    template <typename Derived>
    struct ManagerItemTpl : public internal::CRTP<Derived>
    {
        using MetaManager_t = typename traits<Derived>::MetaManager_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        ManagerItemTpl(const std::string &name_, const Model_t &model_, const bool active_ = true)
            : name(name_), model(model_), active(active_)
        {
        }

        std::string name;
        Model_t model;
        bool active;

    }; // struct ManagerItemTpl

    template <typename Derived>
    class ManagerDataBase : public internal::CRTP<Derived>
    {
    public:
        using MetaManager_t = typename traits<Derived>::MetaManager_t;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        DataContainer_t items;

    protected:
        template <typename DataCollector>
        inline ManagerDataBase(const ModelManager_t &model_manager, MemoryArena &arena, DataCollector *const collector)
        {
            items.clear();
            for (typename ModelContainer_t::const_iterator it = model_manager.get_items().begin();
                 it != model_manager.get_items().end();
                 ++it)
            {
                const Item_t &item = it->second;
                items.insert(std::make_pair(item.name, item.model.createData(arena, collector)));
            }
        }

        inline ManagerDataBase(const ManagerDataBase &clone) : items(clone.items) {}

        inline ManagerDataBase &operator=(const ManagerDataBase &clone)
        {
            items = clone.items;
            return *this;
        }

    }; // class ManagerDataBase

    template <typename Derived>
    class ManagerModelBase : public internal::CRTP<Derived>
    {
    public:
        using MetaManager_t = typename traits<Derived>::MetaManager_t;
        using Collection_t = typename traits<MetaManager_t>::Collection_t;
        using ModelManager_t = typename traits<MetaManager_t>::ModelManager_t;
        using DataManager_t = typename traits<MetaManager_t>::DataManager_t;

        using Meta_t = typename traits<MetaManager_t>::Meta_t;
        using Model_t = typename traits<MetaManager_t>::Model_t;
        using Data_t = typename traits<MetaManager_t>::Data_t;

        using Item_t = typename traits<MetaManager_t>::Item_t;

        using ModelContainer_t = typename traits<MetaManager_t>::ModelContainer_t;
        using DataContainer_t = typename traits<MetaManager_t>::DataContainer_t;

        template <typename StateVectorType, typename ControlVectorType>
        void calc(DataManager_t &data,
                  const Eigen::MatrixBase<StateVectorType> &x,
                  const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calc(data, x, u);
        }

        template <typename StateVectorType>
        void calc(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calc(data, x);
        }

        template <typename StateVectorType, typename ControlVectorType>
        void calcDiff(DataManager_t &data,
                      const Eigen::MatrixBase<StateVectorType> &x,
                      const Eigen::MatrixBase<ControlVectorType> &u) const
        {
            this->derived().calcDiff(data, x, u);
        }

        template <typename StateVectorType>
        void calcDiff(DataManager_t &data, const Eigen::MatrixBase<StateVectorType> &x) const
        {
            this->derived().calcDiff(data, x);
        }

        template <typename DataCollector>
        DataManager_t createData(MemoryArena &arena, DataCollector *const collector) const
        {
            return this->derived().createData(arena, collector);
        }

        template <typename... Args>
        void addItem(const std::string &name, const Model_t &model, Args &&...args)
        {
            this->derived().addItemImpl(name, model, std::forward<Args>(args)...);
        }

        template <typename... Args>
        void addItemImpl(const std::string &name, const Model_t &model, Args &&...args)
        {
            // Extract active status from constructed item.
            // Allows us to keep the default active status but also support 'overloaded' addItemImpl methods that may
            // have different Item constructor signatures
            auto item = Item_t(name, model, std::forward<Args>(args)...);
            const bool active = item.active;

            std::pair<typename ModelContainer_t::iterator, bool> ret =
                items_.insert(std::make_pair(name, std::move(item)));
            if (ret.second == false)
            {
                std::cerr << "Warning: we couldn't add the " << name << " item, it already existed." << std::endl;
            }
            else if (active)
            {
                const int model_n = get_model_n(model);
                active_dim_ += model_n;
                total_dim_ += model_n;
                active_set_.insert(name);
            }
            else if (!active)
            {
                total_dim_ += get_model_n(model);
                inactive_set_.insert(name);
            }
        }

        void removeItem(const std::string &name) { this->derived().removeItemImpl(name); }

        void removeItemImpl(const std::string &name)
        {
            typename ModelContainer_t::iterator it = items_.find(name);
            if (it != items_.end())
            {
                active_dim_ -= get_model_n(it->second.model);
                total_dim_ -= get_model_n(it->second.model);
                items_.erase(it);
                inactive_set_.erase(name);
            }
            else
            {
                std::cerr << "Warning: we couldn't remove the " << name << " item, it doesn't exist." << std::endl;
            }
        }

        void changeItemStatus(const std::string &name, bool active)
        {
            this->derived().changeItemStatusImpl(name, active);
        }

        void changeItemStatusImpl(const std::string &name, bool active)
        {
            typename ModelContainer_t::iterator it = items_.find(name);
            if (it != items_.end())
            {
                if (active && !it->second.active)
                {
                    active_dim_ += get_model_n(it->second.model);
                    active_set_.insert(name);
                    inactive_set_.erase(name);
                    it->second.active = active;
                }
                else if (!active && it->second.active)
                {
                    active_dim_ -= get_model_n(it->second.model);
                    active_set_.erase(name);
                    inactive_set_.insert(name);
                    it->second.active = active;
                }
                it->second.active = active;
            }
            else
            {
                std::cerr << "Warning: we couldn't change the status of the " << name << " item, it doesn't exist."
                          << std::endl;
            }
        }

        const ModelContainer_t &get_items() const { return items_; }
        const std::set<std::string> &get_active_set() const { return active_set_; }
        const std::set<std::string> &get_inactive_set() const { return inactive_set_; }

        bool get_item_status(const std::string &name) const
        {
            typename ModelContainer_t::const_iterator it = items_.find(name);
            if (it != items_.end())
            {
                return it->second.active;
            }
            else
            {
                std::cerr << "Warning: we couldn't get the status of the " << name << " item, it doesn't exist."
                          << std::endl;
                return false;
            }
        }

        const DimensionTpl<Eigen::Dynamic> &get_n_active_dim() const { return active_dim_; }
        int get_n_active() const { return active_dim_.value(); }
        const DimensionTpl<Eigen::Dynamic> &get_n_total_dim() const { return total_dim_; }
        int get_n_total() const { return total_dim_.value(); }

        int get_model_n(const Model_t &model) const { return this->derived().get_model_n(model); }

    protected:
        inline ManagerModelBase()
        {
            items_.clear();
            active_set_.clear();
            inactive_set_.clear();
            active_dim_.set_value(0);
            total_dim_.set_value(0);
        }

        inline ManagerModelBase(const ManagerModelBase &clone)
            : items_(clone.items_),
              active_set_(clone.active_set_),
              inactive_set_(clone.inactive_set_),
              active_dim_(clone.active_dim_),
              total_dim_(clone.total_dim_)
        {
        }

        inline ManagerModelBase &operator=(const ManagerModelBase &clone)
        {
            items_ = clone.items_;
            active_set_ = clone.active_set_;
            inactive_set_ = clone.inactive_set_;
            active_dim_ = clone.active_dim_;
            total_dim_ = clone.total_dim_;
            return *this;
        }

        ModelContainer_t items_;

        std::set<std::string> active_set_;
        std::set<std::string> inactive_set_;

        DimensionTpl<Eigen::Dynamic> active_dim_;
        DimensionTpl<Eigen::Dynamic> total_dim_;

    }; // class ManagerModelBase

} // namespace galileo

#endif // __galileo_common_container_manager_base_hpp__
