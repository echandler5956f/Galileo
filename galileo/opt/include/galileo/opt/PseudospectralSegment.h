#pragma once

#include "galileo/opt/NLPData.h"
#include "galileo/opt/States.h"
#include "galileo/opt/ProblemData.h"
#include "galileo/opt/LagrangePolynomial.h"
#include "galileo/tools/CasadiTools.h"
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
             * @brief Construct a new Pseudospectral Segment object.
             *
             * @param problem Pointer to the problem data
             * @param F Dynamics function
             * @param L Cost function
             * @param st_m Pointer to the state indices helper
             * @param d Polynomial degree
             * @param optimize_dt Flag to optimize the time step
             * @param knot_num Number of knots in the segment
             *
             */
            PseudospectralSegment(std::shared_ptr<GeneralProblemData> problem, casadi::Function F, casadi::Function L, std::shared_ptr<States> st_m, int d, bool optimize_dt, int knot_num);

            /**
             * @brief Build the function graph.
             *
             * @param G Vector of constraint data
             * @param Wdata Decision bound and initial guess data for the state and input
             */
            void InitializeExpressionGraph(std::vector<ConstraintData> G, std::shared_ptr<DecisionData> Wdata, casadi::Dict casadi_opts);

            /**
             * @brief Create all the knot segments.
             *
             * @param X0_sym Global starting state to integrate from
             * @param Xf_sym Global final state
             *
             */
            void InitializeKnotSegments(casadi::MX X0_sym, casadi::MX Xf_sym);

            /**
             * @brief Update a pseudospectral segment.
             *
             * @param segment_offset The starting time of the segment w.r.t the global time
             * @param h The time step
             * @param X0_param The state to deviate from
             * @param Xf_param The final state
             * @param generate_guess Flag to generate an initial guess
             */
            void Update(const double segment_offset, const double h, casadi::DM X0_param, casadi::DM Xf_param, bool generate_guess = true);

            /**
             * @brief Evaliate the expressions with the actual decision variables.
             *
             */
            void EvaluateExpressionGraph();

            /**
             * @brief Extract the solution from the decision variable vector.
             *
             * @param w Decision variable vector
             * @param p0 Initial guess for the parameters
             *
             * @return casadi::DMVector Solution values
             */
            casadi::DMVector ExtractSolution(casadi::DM &w, casadi::DM &p0) const;

            /**
             * @brief Updates the NLP data with the bounds, initial guesses, and parameters values.
             *
             * This function takes in an NLPInputData object, and fills its members with the members of this->nlp_in_data_.
             *
             * @param nlp_in_data The data to be filled
             * @param update_guess Flag to update the initial guess using the initial guess function.
             * If false, the initial guess must be filled elsewhere
             */
            void UpdateNLPInputData(NLPInputData &nlp_in_data, bool update_guess = true);

            /**
             * @brief Fills the NLP data with the decision variables, constraints, and costs.
             *
             * This function takes in an NLPProblemData object, and fills its members with the members of this->nlp_prob_data_.
             *
             * @param nlp_prob_data The data to be filled
             */
            void FillNLPProblemData(NLPProblemData &nlp_prob_data);

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
             * @brief Get the state times vectors.
             *
             * @return PseudospectralTimes The state time vectors
             */
            PseudospectralTimes getStateTimes() const
            {
                return dXtimes_;
            }

            /**
             * @brief Get the input time vectors.
             *
             * @return PseudospectralTimes The input time vectors
             */
            PseudospectralTimes getInputTimes() const
            {
                return Utimes_;
            }

            /**
             * @brief Get the state segment times vector.
             *
             * @return casadi::DM The state segment time vector
             */
            casadi::DM getStateSegmentTimes() const
            {
                return dXtimes_.segment_times;
            }

            /**
             * @brief Get the input segment time vector.
             *
             * @return casadi::DM The input segment time vector
             */
            casadi::DM getInputSegmentTimes() const
            {
                return Utimes_.segment_times;
            }

            /**
             * @brief Get the times in the order that the decision variables are arranged for this segment.
             *
             * @return casadi::DM The stacked decision variable times
             */
            casadi::DMVector getSegmentDecisionVariableTimes() const
            {
                return casadi::DMVector{vertcat(casadi::DMVector{dXtimes_.knot_times, dXtimes_.collocation_times}), vertcat(casadi::DMVector{Utimes_.knot_times, Utimes_.collocation_times})};
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

            /**
             * @brief Get the starting and ending index of the parameters in p corresponding to this segment.
             *
             * @return tuple_size_t Range of the parameters in p
             */
            tuple_size_t getRangeParameters() const
            {
                return p_range_;
            }

        private:
            /**
             * @brief Initialize the relevant expressions.
             *
             * @param d Polynomial degree
             */
            void InitializeExpressionVariables(int d);

            /**
             * @brief Initialize the vector of segment times which constraints are evaluated at.
             *
             * @param segment_offset The starting time of the segment w.r.t the global time
             * @param h The time step
             */
            void UpdateTimeVectors(const double segment_offset, const double h);

            /**
             * @brief Update the bounds with the current time segment.
             *
             */
            void UpdateBounds();

            /**
             * @brief Reset w0, p0, lbw, ubw, lbg, and ubg.
             *
             */
            void ResetNLPInputData();

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
             * @brief Flag to optimize the time step. Not implemented yet.
             * 
             */
            bool optimize_dt_ = false;

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

                /**
                 * @brief Time step used to evaluate the expression graphs.
                 *
                 */
                casadi::MX dt;

                /**
                 * @brief Symbolic global initial state.
                 *
                 */
                casadi::MX X0_sym_;

                /**
                 * @brief
                 *
                 */
                casadi::MX Xf_sym_;
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

                /**
                 * @brief Time step used to build the expression graphs.
                 *
                 */
                casadi::SX dt;

                /**
                 * @brief
                 *
                 */
                casadi::SX Xf_global;
            };

            ExpressionVariables expr_v_;

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
             * @brief Stores the symbolic decision variables, constraints, and costs for this segment for the NLP problem.
             *
             */
            NLPProblemData nlp_prob_data_;

            /**
             * @brief Stores the numeric initial guesses, bounds, and parameters for this segment for the NLP problem.
             *
             */
            NLPInputData nlp_in_data_;

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

                std::vector<casadi::Function> lbg_maps;
                std::vector<casadi::Function> ubg_maps;

                casadi::Function w0_func;
                casadi::Function lbw_func;
                casadi::Function ubw_func;

                std::vector<tuple_size_t> ranges_G;
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
             * @brief Starting and ending index of the decision variables in w corresponding to this segment.
             *
             */
            tuple_size_t w_range_;

            /**
             * @brief Starting and ending index of the parameters in p corresponding to this segment.
             *
             */
            tuple_size_t p_range_;
        };
    }
}