#pragma once

#include "galileo/opt/Segment.h"
#include "galileo/opt/LagrangePolynomial.h"

namespace galileo
{
    namespace opt
    {

        /**
         * @brief PseudospectalSegment class.
         *
         */
        class PseudospectralSegment : public Segment
        {
        public:
            /**
             * @brief Construct a new Pseudospectral Segment object.
             *
             * @param problem Pointer to the problem data
             * @param F Dynamics function
             * @param st_m_ Pointer to the state indices helper
             * @param d Polynomial degree
             * @param knot_num_ Number of knots in the segment
             * @param h_ Period of each knot segment
             *
             */
            PseudospectralSegment(std::shared_ptr<GeneralProblemData> problem, casadi::Function F, std::shared_ptr<States> st_m_, int d, int knot_num_, double h_);

            /**
             * @brief Initialize the relevant expressions.
             *
             * @param d Polynomial degree
             */
            void initializeExpressionVariables(int d);

            /**
             * @brief Initialize the vector of local times which constraints are evaluated at.
             *
             * @param global_times Vector of global times
             */
            void initializeLocalTimeVector(casadi::DM &global_times);

            /**
             * @brief Initialize the vector of times which coincide to the decision variables U occur at.
             *
             */
            void initializeInputTimeVector();

            /**
             * @brief Create all the knot segments.
             *
             * @param x0_global Global starting state to integrate from (used for initial guess)
             * @param x0_local Local starting state to integrate from
             *
             */
            void initializeKnotSegments(casadi::DM x0_global, casadi::MX x0_local);

            /**
             * @brief Build the function graph.
             *
             * @param G Vector of constraint data
             * @param Wx Decision bound and initial guess data for the state
             * @param Wu Decision bound and initial guess data for the input
             */
            void initializeExpressionGraph(std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu);

            /**
             * @brief Evaluate the expressions with the actual decision variables.
             *
             * @param J0 Accumulated cost so far
             * @param w Decision variable vector to fill
             * @param g Constraint vector to fill
             */
            void evaluateExpressionGraph(casadi::MX &J0, casadi::MXVector &w, casadi::MXVector &g);

            /**
             * @brief Extract the solution from the decision variable vector.
             *
             * @param w Decision variable vector
             * @return casadi::MX Solution values
             */
            casadi::MXVector extractSolution(casadi::MX &w) const;

            /**
             * @brief Get the initial state.
             *
             * @return casadi::MX The initial state
             */
            casadi::MX getInitialState() const;

            /**
             * @brief Get the initial state deviant.
             *
             * @return casadi::MX The initial state deviant
             */
            casadi::MX getInitialStateDeviant() const;

            /**
             * @brief Get the final state deviant.
             *
             * @return casadi::MX The final state deviant
             */
            casadi::MX getFinalStateDeviant() const;

            /**
             * @brief Get the actual final state.
             *
             * @return casadi::MX The final state.
             */
            casadi::MX getFinalState() const;

            /**
             * @brief Get the local times vector.
             *
             * @return casadi::DM The local times vector
             */
            casadi::DM getLocalTimes() const;

            /**
             * @brief Fills the lower bounds on decision variable (lbw) and upper bounds on decision variable (ubw) vectors with values.
             *
             * This function takes in two vectors, lbw and ubw, and fills them with values.
             * The filled values represent the constraint lower and upper bounds on decision variables.
             *
             * @param lbw The vector to be filled with lower bound values on decision variables.
             * @param ubw The vector to be filled with upper bound values on decision variables.
             */
            void fill_lbw_ubw(std::vector<double> &lbw, std::vector<double> &ubw);

            /**
             * @brief Fills the lower bounds on general constraints (lbg) and upper bounds on general constraints (ubg) vectors with values.
             *
             * This function takes in two vectors, lbg and ubg, and fills them with values.
             * The filled values represent the general constraint lower and upper bounds.
             *
             * @param lbg The vector to be filled with general lower bound values.
             * @param ubg The vector to be filled with general upper bound values.
             */
            void fill_lbg_ubg(std::vector<double> &lbg, std::vector<double> &ubg);

            /**
             * @brief Fills the initial guess vector (w0) with values.
             *
             * This function takes in a vector, w0, and fills it with values.
             * The filled values represent the initial guess for the decision variables.
             *
             * @param w0 The vector to be filled with initial guess values.
             */
            void fill_w0(std::vector<double> &w0) const;

            /**
             * @brief Returns the starting and ending index in w.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_decision_variables() const;

            /**
             * @brief Returns the starting and ending index in g (call after evaluate_expression_graph!).
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_constraint_expressions() const;

            /**
             * @brief Returns the starting and ending index in lbg/ubg.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_constraint_bounds() const;

            /**
             * @brief Returns the starting and ending index in lbw/ubw.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_decision_bounds() const;

            const std::shared_ptr<LagrangePolynomial> get_dXPoly() const
            {
                return std::make_shared<LagrangePolynomial>(dX_poly);
            }

            const std::shared_ptr<LagrangePolynomial> get_UPoly() const
            {
                return std::make_shared<LagrangePolynomial>(U_poly);
            }

            int getDegree() const
            {
                return dX_poly.d;
            }

            /**
             * @brief Get the knot num
             *
             * @return int
             */
            int getKnotNum() const
            {
                return knot_num;
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
            casadi::MX processVector(casadi::MXVector &vec) const;

            /**
             * @brief Process the offset vector by removing the first element and concatenating the remaining elements horizontally.
             *
             * @param vec The input offset vector.
             * @return The processed offset vector.
             */
            casadi::MX processOffsetVector(casadi::MXVector &vec) const;

            /**
             * @brief Input polynomial. Helper object to store polynomial information for the input.
             *
             */
            LagrangePolynomial U_poly;

            /**
             * @brief State polynomial. Helper object to store polynomial information for the state.
             *
             */
            LagrangePolynomial dX_poly;

            /**
             * @brief Collocation state decision variables.
             *
             */
            casadi::MXVector dXc_var_vec;

            /**
             * @brief Collocation input decision variables.
             *
             */
            casadi::MXVector U_var_vec;

            /**
             * @brief Collocation input decision expressions at the state collocation points.
             *
             * Decision variables of control and state are potentially approximated by different degree polynomials.
             *
             */
            casadi::SXVector U_at_c_vec;

            /**
             * @brief Collocation state decision expressions at the collocation points.
             *
             */
            casadi::SXVector x_at_c_vec;

            /**
             * @brief Knot point deviants state decision variables.
             *
             */
            casadi::MXVector dX0_var_vec;

            /**
             * @brief Knot point state expressions (integral functions of the deviants).
             *
             */
            casadi::MXVector X0_var_vec;

            /**
             * @brief casadi::Function map for extracting the solution from the ocp solution vector.
             *
             */
            casadi::Function sol_map_func;

            /**
             * @brief Solution function.
             *
             */
            casadi::Function get_sol_func;

            /**
             * @brief Implicit discrete-time function map. This function map returns the vector of collocation equations
                necessary to match the derivative defect between the approximated dynamics and actual system
                dynamics.
             *
             */
            casadi::Function collocation_constraint_map;

            /**
             * @brief Implicit discrete-time function map. The map which matches the approximated final state expression with the initial
                state of the next segment.
             *
             */
            casadi::Function xf_constraint_map;

            /**
             * @brief Implicit discrete-time function map. The accumulated cost across all the knot segments found using quadrature rules.
             *
             */
            casadi::Function q_cost_fold;

            /**
             * @brief User defined constraints, which are functions with certain bounds associated with them.
             *
             */
            std::vector<casadi::Function> general_constraint_maps;

            /**
             * @brief Lower bounds associated with the general constraint maps.
             *
             */
            casadi::DM general_lbg;

            /**
             * @brief Upper bounds associated with the general constraint maps.
             *
             */
            casadi::DM general_ubg;

            /**
             * @brief Lower bounds associated with the decision variable constraint maps.
             *
             */
            casadi::DM general_lbw;

            /**
             * @brief Upper bounds associated with the decision variable constraint maps.
             *
             */
            casadi::DM general_ubw;

            /**
             * @brief Initial guess for associated with this segment
             */
            casadi::DM w0;

            /**
             * @brief Integrator function.
             *
             */
            casadi::Function Fint;

            /**
             * @brief Difference function.
             *
             */
            casadi::Function Fdif;

            /**
             * @brief Dynamics function.
             *
             */
            casadi::Function F;

            /**
             * @brief Cost function.
             *
             */
            casadi::Function L;

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
             * @brief Accumulator expression used to build the expression graphs.
             *
             */
            casadi::SX Lc;

            /**
             * @brief Helper for indexing the state variables.
             *
             */
            std::shared_ptr<States> st_m;

            /**
             * @brief Number of knot segments.
             *
             */
            int knot_num;

            /**
             * @brief Vector of local times including knot points. Note that this coincides with the times for the decision variables of x.
             *
             */
            casadi::DM local_times;

            /**
             * @brief Vector of unique times for this segment, including terminal but not including initial.
             * Used for bound constraint evaluation, so that we don't overconstrain variables which occur at the same time
             * e.g, x0 and xf of adjacent segments
             *
             * In the frame of the global times (NOT the local times)
             *
             */
            casadi::DM collocation_times;

            /**
             * @brief Vector of times for the knot points for this segment.
             *
             * In the frame of the global times (NOT the local times)
             */
            casadi::DM knot_times;

            /**
             * @brief Vector of times for the decision variables of u.
             */
            casadi::DM u_times;
        };
    }
}