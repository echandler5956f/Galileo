#pragma once

#include "galileo/variables/States.h"
#include "galileo/variables/Constraint.h"
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <cassert>
#include <bits/stdc++.h>
#include <memory>

using namespace casadi;

namespace galileo
{
    namespace variables
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
            const SX lagrange_interpolation(double t, const SXVector terms);

            /**
             * @brief Degree of the polynomial.
             *
             */
            int d;

            /**
             * @brief The roots of the polynomial.
             *
             */
            Eigen::VectorXd tau_root;

            /**
             * @brief Quadrature coefficients.
             *
             */
            Eigen::VectorXd B;

            /**
             * @brief Collocation coeffficients.
             *
             */
            Eigen::MatrixXd C;

            /**
             * @brief Continuity coefficients.
             *
             */
            Eigen::VectorXd D;
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
             * @param st_m_ Pointer to the state indices helper
             * @param Fint_ Integrator function
             */
            PseudospectralSegment(int d, int knot_num_, double h_, std::shared_ptr<States> st_m_, Function &Fint_);

            /**
             * @brief Initialize the relevant expressions.
             *
             * @param d Polynomial degree
             */
            void initialize_expression_variables(int d);

            /**
             * @brief Initialize the vector of all times which constraints are evaluated at.
             *
             */
            void initialize_time_vector();

            /**
             * @brief Fill all times with the time vector from this segment.
             *
             * @param all_times
             */
            void fill_times(std::vector<double> &all_times);

            /**
             * @brief Create all the knot segments.
             *
             * @param x0 Starting state to integrate from. Can be a constant
             *
             */
            void initialize_knot_segments(SX x0);

            /**
             * @brief Build the function graph.
             *
             * @param F Function for the system dynamics
             * @param L Integrated cost
             * @param G Vector of constraint data
             */
            void initialize_expression_graph(Function &F, Function &L, std::vector<std::shared_ptr<ConstraintData>> G);

            /**
             * @brief Evaluate the expressions with the actual decision variables.
             *
             * @param J0 Accumulated cost so far
             * @param g Constraint vector to fill
             */
            void evaluate_expression_graph(SX &J0, SXVector &g);

            /**
             * @brief Get the initial state.
             *
             * @return SX The initial state
             */
            SX get_initial_state();

            /**
             * @brief Get the initial state deviant.
             *
             * @return SX The initial state deviant
             */
            SX get_initial_state_deviant();

            /**
             * @brief Get the final state deviant.
             *
             * @return SX The final state deviant
             */
            SX get_final_state_deviant();

            /**
             * @brief Get the actual final state.
             *
             * @return SX The final state.
             */
            SX get_final_state();

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
             * @brief Fills the lower bounds on decision variable (lbx) and upper bounds on decision variable (ubx) vectors with values.
             *
             * This function takes in two vectors, lbx and ubx, and fills them with values.
             * The filled values represent the constraint lower and upper bounds on decision variables.
             *
             * @param lbx The vector to be filled with lower bound values on decision variables.
             * @param ubx The vector to be filled with upper bound values on decision variables.
             */
            void fill_lbx_ubx(std::vector<double> &lbx, std::vector<double> &ubx);

            /**
             * @brief Fills the given SXVector with values.
             *
             * This function fills the provided SXVector with values.
             *
             * @param w The SXVector to be filled.
             */
            void fill_w(SXVector &w);

            /**
             * @brief Returns the starting and ending index in w (call after fill_w!).
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_decision_variables();

            /**
             * @brief Returns the starting and ending index in g (call after evaluate_expression_graph!).
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_constraint_expressions();

            /**
             * @brief Returns the starting and ending index in bg (call after fill_lbg_ubg!). This should match get_range_idx_constraint_expressions.
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_bg();

            /**
             * @brief Returns the starting and ending index in bx (call after fill_lbx_ubx!).
             *
             * @return tuple_size_t The range of indices
             */
            tuple_size_t get_range_idx_bx();

        private:
            /**
             * @brief Collocation state decision variables.
             *
             */
            SXVector dXc_var_vec;

            /**
             * @brief Collocation input decision variables.
             *
             */
            SXVector U_var_vec;

            /**
             * @brief Collocation input decision expressions at the state collocation points.
             *
             * Decision variables of control and state are potentially approximated by different degree polynomials.
             *
             */
            SXVector U_at_c_vec;

            /**
             * @brief Knot point deviants state decision variables.
             *
             */
            SXVector dX0_var_vec;

            /**
             * @brief Knot point state expressions (integral functions of the deviants).
             *
             */
            SXVector X0_var_vec;

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
            DM general_lbx;

            /**
             * @brief Upper bounds associated with the decision variable constraint maps.
             *
             */
            DM general_ubx;

            /**
             * @brief Integrator function.
             *
             */
            Function Fint;

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
             * @brief Vector of all times.
             *
             */
            DM times;

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
             * @brief Starting and ending index of the constraint expressions in g corresponding to this segment. This should match lbg_ubg_range.
             *
             */
            tuple_size_t g_range;

            /**
             * @brief Starting and ending index of the constraint expressions in bx corresponding to this segment. This should match lbx_ubx_range.
             *
             */
            tuple_size_t bx_range;

            /**
             * @brief Starting and ending index of the bounds in lbg/ubg corresponding to this segment. This should match g_range.
             *
             */
            tuple_size_t lbg_ubg_range;

            /**
             * @brief Starting and ending index of the bounds in lbx/ubx corresponding to this segment. This should match bx_range.
             *
             */
            tuple_size_t lbx_ubx_range;
        };
    }
}