#ifndef __galileo_core_constraint_collections_hpp__
#define __galileo_core_constraint_collections_hpp__

#include "galileo/core/constraint/fwd.hpp"

#include <bitset>
#include <tuple>
#include <utility>
#include <type_traits>

namespace galileo
{

    // -----------------------------------------------------------------------------
    // 5. The ConstraintCollection class
    //    - Stores constraints in a tuple
    //    - Maintains a bitset to track which constraints are active
    //    - Provides type-based access (activate<T>(), deactivate<T>(), get<T>(), etc.)
    //    - Uses perfect forwarding in constructor
    //    - Demonstrates a "calcAllActive" that zips constraints and data
    // -----------------------------------------------------------------------------

    template <typename... Constraints>
    class ConstraintCollection
    {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        // By default, either all active or all inactive.
        // Here, let's set them all active by default:
        ConstraintCollection()
        {
            activeMask_.set(); // all active
        }

        // Perfect forwarding constructor:
        // allow the user to provide initial constraints as constructor arguments
        // e.g., ConstraintCollection(SomeSpecificConstraint<double>{}, AnotherConstraint<double>{})
        template <typename... CtorConstraints,
                  typename = std::enable_if_t<(sizeof...(CtorConstraints) == num_constraints_)>>
        explicit ConstraintCollection(CtorConstraints &&...cs)
            : constraints_(std::forward<CtorConstraints>(cs)...)
        {
            activeMask_.set(); // default all active
        }

        // ------------------------------------------------------------------------
        // Activation / Deactivation by type
        // ------------------------------------------------------------------------

        template <typename T>
        void activate()
        {
            constexpr std::size_t idx = index_of<T, Constraints...>::value;
            activeMask_.set(idx, true);
        }

        template <typename T>
        void deactivate()
        {
            constexpr std::size_t idx = index_of<T, Constraints...>::value;
            activeMask_.set(idx, false);
        }

        template <typename T>
        bool isActive() const
        {
            constexpr std::size_t idx = index_of<T, Constraints...>::value;
            return activeMask_.test(idx);
        }

        // ------------------------------------------------------------------------
        // Direct, type-based access to constraint objects:
        //    collection.get<SomeSpecificConstraint<double>>().setBounds(...)
        // ------------------------------------------------------------------------
        template <typename T>
        T &get()
        {
            return std::get<T>(constraints_);
        }
        template <typename T>
        const T &get() const
        {
            return std::get<T>(constraints_);
        }

        // ------------------------------------------------------------------------
        // calcAllActive: calls "calc" on each constraint that is active,
        // pairing with the corresponding data object in dataTuple.
        //
        // dataTuple must have the same length (# of elements) as constraints_.
        // For example:
        //   std::tuple<SomeSpecificConstraintData, AnotherConstraintData>
        // ------------------------------------------------------------------------
        template <typename DataTuple, typename State, typename Control>
        void calcAllActive(DataTuple &dataTuple,
                           const State &xs,
                           const Control &us)
        {
            static_assert(std::tuple_size<DataTuple>::value == sizeof...(Constraints),
                          "Data tuple must have same number of elements as constraints.");

            calcAllActiveImpl(dataTuple, xs, us, std::make_index_sequence<sizeof...(Constraints)>{});
        }

    private:
        // Implementation detail: expand over indices [0..num_constraints_-1]
        template <typename DataTuple, typename State, typename Control, std::size_t... Is>
        void calcAllActiveImpl(DataTuple &dataTuple,
                               const State &xs,
                               const Control &us,
                               std::index_sequence<Is...>)
        {
            // We'll fold over doOneActive<Is> calls
            (doOneActive<Is>(dataTuple, xs, us), ...);
        }

        template <std::size_t I, typename DataTuple, typename State, typename Control>
        void doOneActive(DataTuple &dataTuple,
                         const State &xs,
                         const Control &us)
        {
            // If constraint I is active, call calc
            if (activeMask_.test(I))
            {
                auto &constraint = std::get<I>(constraints_);
                auto &data = std::get<I>(dataTuple);
                // Now call the CRTP-based calc
                constraint.calc(data, xs, us);
            }
        }

        std::tuple<Constraints...> constraints_;   // The constraints
        std::bitset<sizeof...(Constraints)> activeMask_; // Which constraints are active?
    };

} // namespace galileo

#endif // __galileo_core_constraint_collections_hpp__
