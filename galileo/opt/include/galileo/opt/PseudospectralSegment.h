#pragma once

#include "galileo/opt/States.h"
#include "galileo/opt/Constraint.h"
#include <vector>
#include <string>
#include <cassert>
#include <memory>

using namespace casadi;

namespace galileo
{
    namespace opt
    {

        using tuple_size_t = std::tuple<std::size_t, std::size_t>;

        /**
         * @brief Helper class for storing polynomial information.
         *
         */
        class LagrangePolynomial
        {
        public:
            /**
             * @brief Construct a new Lagrange Polynomial object.
             *
             */
            LagrangePolynomial(){};

            /**
             * @brief Construct a new Lagrange Polynomial object. Compute and store the coefficients for a given degree and collocation scheme.
             *
             * @param d_ Degree of the polynomial
             * @param scheme Collocation scheme: "radau" or "legendre"
             */
            LagrangePolynomial(int d_, const std::string &scheme = "radau");

            /**
             * @brief Perform symbolic Lagrange Interpolation, which, given a time from the Lagrange time scale, interpolates terms to find the value at time t.
             *
             * @param t Time to interpolate at
             * @param terms Terms at knot points to use for interpolation
             * @return const SX Resultant expression for the symbolic interpolated value
             */
            SX lagrange_interpolation(double t, const SXVector terms) const;

            /**
             * @brief Degree of the polynomial.
             *
             */
            int d;

            /**
             * @brief The roots of the polynomial.
             *
             */
            std::vector<double> tau_root;

            /**
             * @brief Quadrature coefficients.
             *
             */
            std::vector<double> B;

            /**
             * @brief Collocation coeffficients.
             *
             */
            std::vector<std::vector<double>> C;

            /**
             * @brief Continuity coefficients.
             *
             */
            std::vector<double> D;
        };

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
             */
            PseudospectralSegment(){};

            /**
             * @brief Construct a new Pseudospectral Segment object.
             *
             * @param d Polynomial degree
             * @param knot_num_ Number of knots in the segment
             * @param h_ Period of each knot segment
             * @param global_times_ Vector of glbal times
             * @param st_m_ Pointer to the state indices helper
             * @param Fint_ Integrator function
             * @param Fdif_ Difference function
             */
            PseudospectralSegment(int d, int knot_num_, double h_, std::shared_ptr<casadi::DM> global_times_, std::shared_ptr<States> st_m_, Function &Fint_, Function &Fdif_);

            /**
             * @brief Initialize the relevant expressions.
             *
             * @param d Polynomial degree
             */
            void initialize_expression_variables(int d);

            /**
             * @brief Initialize the vector of local times which constraints are evaluated at.
             *
             */
            void initialize_local_time_vector();

            /**
             * @brief Initialize the vector of times which coincide to the decision variables U occur at.
             *
             */
            void initialize_u_time_vector();

            /**
             * @brief Create all the knot segments.
             *
             * @param x0_global Global starting state to integrate from (used for initial guess)
             * @param x0_local Local starting state to integrate from
             *
             */
            void initialize_knot_segments(DM x0_global, MX x0_local);

            /**
             * @brief Build the function graph.
             *
             * @param F Function for the system dynamics
             * @param L Integrated cost
             * @param G Vector of constraint data
             * @param Wx Decision bound and initial guess data for the state
             * @param Wu Decision bound and initial guess data for the input
             */
            void initialize_expression_graph(Function &F, Function &L, std::vector<std::shared_ptr<ConstraintData>> G, std::shared_ptr<DecisionData> Wx, std::shared_ptr<DecisionData> Wu);

            /**
             * @brief Evaluate the expressions with the actual decision variables.
             *
             * @param J0 Accumulated cost so far
             * @param w Decision variable vector to fill
             * @param g Constraint vector to fill
             */
            void evaluate_expression_graph(MX &J0, MXVector &w, MXVector &g);

            /**
             * @brief Extract the solution from the decision variable vector.
             *
             * @param w Decision variable vector
             * @return MX Solution values
             */
            MXVector extract_solution(MX &w) const;

            /**
             * @brief Get the initial state.
             *
             * @return MX The initial state
             */
            MX get_initial_state() const;

            /**
             * @brief Get the initial state deviant.
             *
             * @return MX The initial state deviant
             */
            MX get_initial_state_deviant() const;

            /**
             * @brief Get the final state deviant.
             *
             * @return MX The final state deviant
             */
            MX get_final_state_deviant() const;

            /**
             * @brief Get the actual final state.
             *
             * @return MX The final state.
             */
            MX get_final_state() const;

            /**
             * @brief Get the global times vector.
             *
             * @return std::shared_ptr<casadi::DM> The global times vector
             */
            std::shared_ptr<casadi::DM> get_global_times() const;

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

        private:
            /**
             * @brief Helper function to process a vector of type MX.
             *
             * This function takes a vector of type MX and performs some processing on it.
             * It creates a temporary vector by copying the input vector and removing the last element.
             * The modified vector is then returned.
             *
             * @param vec The input vector of type MX.
             * @return The processed vector of type MX.
             */
            MX processVector(MXVector &vec) const;

            /**
             * @brief Process the offset vector by removing the first element and concatenating the remaining elements horizontally.
             *
             * @param vec The input offset vector.
             * @return The processed offset vector.
             */
            MX processOffsetVector(MXVector &vec) const;

            /**
             * @brief Local initial state.
             *
             */
            MX x0_local;

            /**
             * @brief Global initial state.
             *
             */
            DM x0_global;

            /**
             * @brief Collocation state decision variables.
             *
             */
            MXVector dXc_var_vec;

            /**
             * @brief Collocation input decision variables.
             *
             */
            MXVector U_var_vec;

            /**
             * @brief Collocation input decision expressions at the state collocation points.
             *
             * Decision variables of control and state are potentially approximated by different degree polynomials.
             *
             */
            SXVector U_at_c_vec;

            /**
             * @brief Collocation state decision expressions at the collocation points.
             *
             */
            SXVector x_at_c_vec;

            /**
             * @brief Knot point deviants state decision variables.
             *
             */
            MXVector dX0_var_vec;

            /**
             * @brief Knot point state expressions (integral functions of the deviants).
             *
             */
            MXVector X0_var_vec;

            /**
             * @brief Function map for converting solution to plottable results.
             *
             */
            Function plot_map_func;

            /**
             * @brief Solution function.
             *
             */
            Function get_sol_func;

            /**
             * @brief Implicit discrete-time function map. This function map returns the vector of collocation equations
                necessary to match the derivative defect between the approximated dynamics and actual system
                dynamics.
             *
             */
            Function collocation_constraint_map;

            /**
             * @brief Implicit discrete-time function map. The map which matches the approximated final state expression with the initial
                state of the next segment.
             *
             */
            Function xf_constraint_map;

            /**
             * @brief Implicit discrete-time function map. The accumulated cost across all the knot segments found using quadrature rules.
             *
             */
            Function q_cost_fold;

            /**
             * @brief User defined constraints, which are functions with certain bounds associated with them.
             *
             */
            std::vector<Function> general_constraint_maps;

            /**
             * @brief Lower bounds associated with the general constraint maps.
             *
             */
            DM general_lbg;

            /**
             * @brief Upper bounds associated with the general constraint maps.
             *
             */
            DM general_ubg;

            /**
             * @brief Lower bounds associated with the decision variable constraint maps.
             *
             */
            DM general_lbw;

            /**
             * @brief Upper bounds associated with the decision variable constraint maps.
             *
             */
            DM general_ubw;

            /**
             * @brief Initial guess for associated with this segment
             */
            DM w0;

            /**
             * @brief Integrator function.
             *
             */
            Function Fint;

            /**
             * @brief Difference function.
             *
             */
            Function Fdif;

            /**
             * @brief Collocation states used to build the expression graphs.
             *
             */
            SXVector dXc;

            /**
             * @brief Collocation inputs used to build the expression graphs.
             *
             */
            SXVector Uc;

            /**
             * @brief Knot states deviants used to build the expression graphs.
             *
             */
            SX dX0;

            /**
             * @brief Knot states used to build the expression graphs.
             *
             */
            SX X0;

            /**
             * @brief Accumulator expression used to build the expression graphs.
             *
             */
            SX Lc;

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
             * @brief Vector of global times including knot points.
             *
             */
            std::shared_ptr<DM> global_times;

            /**
             * @brief Vector of local times including knot points. Note that this coincides with the times for the decision variables of x.
             *
             */
            DM local_times;

            /**
             * @brief Vector of unique times for this segment, including terminal but not including initial.
             * Used for bound constraint evaluation, so that we don't overconstrain variables which occur at the same time
             * e.g, x0 and xf of adjacent segments
             *
             * In the frame of the global times (NOT the local times)
             *
             */
            DM collocation_times;

            /**
             * @brief Vector of times for the knot points for this segment.
             *
             * In the frame of the global times (NOT the local times)
             */
            DM knot_times;

            /**
             * @brief Vector of times for the decision variables of u.
             */
            DM u_times;

            /**
             * @brief Period of EACH KNOT SEGMENT within this pseudospectral segment.
             *
             */
            double h;

            /**
             * @brief Total period (helper variable calculated from h and knot_num).
             *
             */
            double T;

            /**
             * @brief Starting and ending index of the decision variables in w corresponding to this segment.
             *
             */
            tuple_size_t w_range;

            /**
             * @brief Starting and ending index of the constraint expressions in g corresponding to this segment.
             *
             */
            tuple_size_t g_range;

            /**
             * @brief Starting and ending index of the bounds in lbg/ubg corresponding to this segment.
             *
             */
            tuple_size_t lbg_ubg_range;

            /**
             * @brief Starting and ending index of the bounds in lbw/ubw corresponding to this segment.
             *
             */
            tuple_size_t lbw_ubw_range;
        };
    }
}