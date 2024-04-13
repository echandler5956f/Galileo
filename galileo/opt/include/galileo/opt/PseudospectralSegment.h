#pragma once

#include "galileo/opt/NLPData.h"
#include "galileo/opt/States.h"
#include "galileo/opt/ProblemData.h"
#include "galileo/opt/LagrangePolynomial.h"
#include "galileo/tools/CasadiConversions.h"
#include <cassert>
#include <chrono>

namespace galileo
{
    namespace opt
    {
        using tuple_size_t = std::tuple<size_t, size_t>;

        /**
         * @brief PseudospectalSegment class.
         *
         */
        class PseudospectralSegment
        {
        public:
            /**
             * @brief Construct a new Pseudospectral Segment object.
             *
             * @param problem Pointer to the problem data
             * @param F Dynamics function
             * @param L Cost function
             * @param st_m Pointer to the state indices helper
             * @param d Polynomial degree
             * @param knot_num Number of knots in the segment
             * @param h Period of each knot segment
             *
             */
            PseudospectralSegment(std::shared_ptr<GeneralProblemData> problem, casadi::Function F, casadi::Function L, std::shared_ptr<States> st_m, int d, int knot_num, double h);

            /**
             * @brief Initialize the relevant expressions.
             *
             * @param d Polynomial degree
             */
            void InitializeExpressionVariables(int d);

            /**
             * @brief Initialize the vector of segment times which constraints are evaluated at.
             *
             * @param trajectory_times Vector of global times
             */
            void InitializeTimeVectors(casadi::DM &trajectory_times);

            /**
             * @brief Create all the knot segments.
             *
             * @param x0_global Global starting state to integrate from (used for initial guess)
             *
             */
            void InitializeKnotSegments(casadi::DM x0_global);

            /**
             * @brief Build the function graph.
             *
             * @param G Vector of constraint data
             * @param Wdata Decision bound and initial guess data for the state and input
             */
            void InitializeExpressionGraph(std::vector<ConstraintData> G, std::shared_ptr<DecisionData> Wdata);

            /**
             * @brief Evaluate the expressions with the actual decision variables.
             *
             */
            void EvaluateExpressionGraph();

            /**
             * @brief Extract the solution from the decision variable vector.
             *
             * @param w Decision variable vector
             * @return casadi::MX Solution values
             */
            casadi::MXVector ExtractSolution(casadi::MX &w) const;

            /**
             * @brief Fills the NLP data with the decision variables, constraints, costs, and bounds.
             *
             * This function takes in an NLPData object, and fills its members with the members of this->nlp_data_.
             *
             * @param lbw The vector to be filled with lower bound values on decision variables.
             */
            void FillNLPData(NLPData &nlp_data);

            /**
             * @brief Get the initial state.
             *
             * @return casadi::MX The initial state
             */
            casadi::MX getInitialState() const
            {
                return ps_vars_.X0_vec.front();
            }

            /**
             * @brief Get the initial state deviant.
             *
             * @return casadi::MX The initial state deviant
             */
            casadi::MX getInitialStateDeviant() const
            {
                return ps_vars_.dX0_vec.front();
            }

            /**
             * @brief Get the actual final state.
             *
             * @return casadi::MX The final state.
             */
            casadi::MX getFinalState() const
            {
                return ps_vars_.X0_vec.back();
            }

            /**
             * @brief Get the final state deviant.
             *
             * @return casadi::MX The final state deviant
             */
            casadi::MX getFinalStateDeviant() const
            {
                return ps_vars_.dX0_vec.back();
            }

            /**
             * @brief Get the segment times vector.
             *
             * @return casadi::DM The segment times vector
             */
            casadi::DM getSegmentTimes() const
            {
                return dXtimes_.segment_times;
            }

            /**
             * @brief Get the knot times vector.
             *
             * @return casadi::DM The knot times vector
             */
            casadi::DM getKnotTimes() const
            {
                return dXtimes_.knot_times;
            }

            /**
             * @brief Get the collocation times vector.
             *
             * @return casadi::DM The collocation times vector
             */
            casadi::DM getCollocationTimes() const
            {
                return dXtimes_.collocation_times;
            }

            /**
             * @brief Get the segment times vector of the input.
             *
             * @return casadi::DM The segment times vector
             */
            casadi::DM getUSegmentTimes() const
            {
                return Utimes_.segment_times;
            }

            /**
             * @brief Get the knot times vector of the input.
             *
             * @return casadi::DM The knot times vector
             */
            casadi::DM getUKnotTimes() const
            {
                return Utimes_.knot_times;
            }

            /**
             * @brief Get the collocation times vector of the input.
             *
             * @return casadi::DM The collocation times vector
             */
            casadi::DM getUCollocationTimes() const
            {
                return Utimes_.collocation_times;
            }

            /**
             * @brief Get the dXPoly object.
             *
             * @return const std::shared_ptr<LagrangePolynomial>
             */
            const std::shared_ptr<LagrangePolynomial> get_dXPoly() const
            {
                return std::make_shared<LagrangePolynomial>(dX_poly_);
            }

            /**
             * @brief Get the UPoly object
             *
             * @return const std::shared_ptr<LagrangePolynomial>
             */
            const std::shared_ptr<LagrangePolynomial> get_UPoly() const
            {
                return std::make_shared<LagrangePolynomial>(U_poly_);
            }

            /**
             * @brief Get the degree
             *
             * @return int The degree
             */
            int getStateDegree() const
            {
                return dX_poly_.d;
            }

            /**
             * @brief Get the input degree
             *
             * @return int The input degree
             */
            int getInputDegree() const
            {
                return U_poly_.d;
            }

            /**
             * @brief Get the knot num
             *
             * @return int The knot num
             */
            int getKnotNum() const
            {
                return knot_num_;
            }

            /**
             * @brief Get the starting and ending index of the decision variables in w corresponding to this segment.
             *
             * @return tuple_size_t Range of the decision variables in w
             */
            tuple_size_t getRangeDecisionVariables() const
            {
                return w_range_;
            }

        private:
            /**
             * @brief Helper function to process a vector of type MX.
             *
             * This function takes a vector of type casadi::MX and performs some processing on it.
             * It creates a temporary vector by copying the input vector and removing the last element.
             * The modified vector is then returned.
             *
             * @param vec The input vector of type casadi::MX.
             * @return The processed vector of type casadi::MX.
             */
            casadi::MX ProcessVector(casadi::MXVector &vec) const;

            /**
             * @brief Process the offset vector by removing the first element and concatenating the remaining elements horizontally.
             *
             * @param vec The input offset vector.
             * @return The processed offset vector.
             */
            casadi::MX ProcessOffsetVector(casadi::MXVector &vec) const;

            /**
             * @brief Input polynomial. Helper object to store polynomial information for the input.
             *
             */
            LagrangePolynomial U_poly_;

            /**
             * @brief State polynomial. Helper object to store polynomial information for the state.
             *
             */
            LagrangePolynomial dX_poly_;

            /**
             * @brief Helper object to store the pseudospectral variables used throughout the segment.
             *
             */
            struct PseudospectralVariables
            {
                /**
                 * @brief Knot point deviants state decision variables.
                 *
                 */
                casadi::MXVector dX0_vec;

                /**
                 * @brief Collocation state decision variables.
                 *
                 */
                casadi::MXVector dXc_vec;

                /**
                 * @brief Knot point input decision variables.
                 *
                 */
                casadi::MXVector U0_vec;

                /**
                 * @brief Collocation input decision variables.
                 *
                 */
                casadi::MXVector Uc_vec;

                /**
                 * @brief Knot point state expressions (integral functions of the deviants).
                 *
                 */
                casadi::MXVector X0_vec;
            };

            PseudospectralVariables ps_vars_;

            /**
             * @brief Helper object to store the expression variables, which are used to build the expression graphs.
             *
             */
            struct ExpressionVariables
            {
                /**
                 * @brief Collocation states used to build the expression graphs.
                 *
                 */
                casadi::SXVector dXc;

                /**
                 * @brief Collocation inputs used to build the expression graphs.
                 *
                 */
                casadi::SXVector Uc;

                /**
                 * @brief Knot states deviants used to build the expression graphs.
                 *
                 */
                casadi::SX dX0;

                /**
                 * @brief Knot states used to build the expression graphs.
                 *
                 */
                casadi::SX X0;

                /**
                 * @brief Knot inputs used to build the expression graphs.
                 *
                 */
                casadi::SX U0;

                /**
                 * @brief Accumulator expression used to build the expression graphs.
                 *
                 */
                casadi::SX Lc;
            };

            ExpressionVariables expr_v_;

            /**
             * @brief Helper object to store the pseudospectral times.
             *
             */
            struct PseudospectralTimes
            {
                casadi::DM segment_times;
                casadi::DM collocation_times;
                casadi::DM knot_times;
            };

            /**
             * @brief Pseudospectral times for the state and state deviants.
             *
             */
            PseudospectralTimes dXtimes_;

            /**
             * @brief Pseudospectral times for the input.
             *
             */
            PseudospectralTimes Utimes_;

            /**
             * @brief Stores the raw decision variables, constraints, costs, and bounds for this segment for the NLP problem.
             *
             */
            NLPData nlp_data_;

            /**
             * @brief casadi::Function map for extracting the solution from the ocp solution vector.
             *
             */
            casadi::Function sol_map_func_;

            /**
             * @brief Solution function.
             *
             */
            casadi::Function get_sol_func_;

            /**
             * @brief Constraint maps and cost folds used to evaluate a pseudospectral segment.
             *
             */
            struct PseudospectralFunctions
            {
                /**
                 * @brief Implicit discrete-time function map. This function map returns the vector of collocation equations
                    necessary to match the derivative defect between the approximated dynamics and actual system
                    dynamics.
                 *
                 */
                casadi::Function col_con_map;

                /**
                 * @brief Implicit discrete-time function map. The map which matches the approximated final state expression with the initial
                    state of the next segment.
                 *
                 */
                casadi::Function xf_con_map;

                /**
                 * @brief Implicit discrete-time function map. The map which matches the approximated final input expression with the initial
                    input of the next segment.
                *
                */
                casadi::Function uf_con_map;

                /**
                 * @brief Implicit discrete-time function map. The accumulated cost across all the knot segments found using quadrature rules.
                 *
                 */
                casadi::Function q_cost_fold;

                /**
                 * @brief User defined constraints, which are functions with certain bounds associated with them.
                 *
                 */
                std::vector<casadi::Function> gcon_maps;
            };

            PseudospectralFunctions ps_funcs_;

            /**
             * @brief Integrator function.
             *
             */
            casadi::Function Fint_;

            /**
             * @brief Difference function.
             *
             */
            casadi::Function Fdiff_;

            /**
             * @brief Dynamics function.
             *
             */
            casadi::Function F_;

            /**
             * @brief Cost function.
             *
             */
            casadi::Function L_;

            /**
             * @brief Helper for indexing the state variables.
             *
             */
            std::shared_ptr<States> st_m_;

            /**
             * @brief Number of knot segments.
             *
             */
            int knot_num_;

            /**
             * @brief Global initial state.
             *
             */
            casadi::DM x0_global_;

            /**
             * @brief Period of EACH KNOT SEGMENT within this pseudospectral segment.
             *
             */
            double h_;

            /**
             * @brief Total period (helper variable calculated from h and knot_num).
             *
             */
            double T_;

            /**
             * @brief Starting and ending index of the time vector corresponding to this segment.
             *
             */
            tuple_size_t t_range_;

            /**
             * @brief Starting and ending index of the decision variables in w corresponding to this segment.
             *
             */
            tuple_size_t w_range_;
        };
    }
}