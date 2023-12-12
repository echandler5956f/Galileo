#include "variables/include/Constraint.h"

namespace acro{
namespace model{
namespace legged{
    /**
     * A builder for the Friction Cone Constraint. 
     * 
     * First Order Approximation -- 
     *      [ 1    0  -mu]
     *      [ 0    1  -mu] 
     *      [-1    0  -mu] * Rotation * GRF  <= 0
     *      [ 0   -1  -mu]
     *      [ 0    0   -1]    
     * 
     *  
     * Second Order Approximation -- 
     *                              [0   0  mu]
     *      LorentzConeConstraint ( [0   1   0] * Rotation * GRF ) 
     *                              [1   0   0]
     * 
     *      or      [0   0  mu] * Rotation * GRF 
     *               - 2norm(   [0   1   0] * Rotation * GRF  )  => 0 
     *                          [1   0   0]
    */
    template <class ProblemData>
    class FrictionConeConstraintBuilder : Constraint<ProblemData>{

        public: 

        enum ApproximationOrder {FIRST_ORDER, SECOND_ORDER};

        FrictionConeConstraintBuilder(double mu, ApproximationOrder approximation_order = ApproximationOrder::FIRST_ORDER) :  
            ConstraintBuilder(), mu_(mu), approximation_order_(approximation_order) {}
        
        /**
         * 
         * @brief Generate flags for each knot point. We set it to all ones, applicable at each knot. 
         * 
         * @param problem_data 
         * @param apply_at 
         */
        void CreateApplyAt(const ProblemData &problem_data, Eigen::VectorXi &apply_at) const override{
            uint num_knots = problem_data.getNumKnots();
            apply_at = Eigen::VectorXi::Constant(num_knots, 1);
        }


        /**
         * @brief Generate bounds for a vector of points. 
         * 
         * For both approximations, each value is less than 0, with no lower bound.
         * 
         * @param problem_data 
         * @param upper_bound 
         * @param lower_bound 
         */
        void CreateBounds(const ProblemData &problem_data, casadi::Function &upper_bound, casadi::Function &lower_bound) const override;

        /**
         * @brief Generate a function to evaluate each point
         * 
         * @param problem_data 
         * @param G
         */
        void CreateFunction(const ProblemData &problem_data, casadi::Function &G) const;

        private:

        void CreateSingleEndEffectorFunction(const std::string& EndEffectorID, const ProblemData &problem_data, casadi::Function &G);
        void CreateSingleEndEffectorBounds(const std::string& EndEffectorID, const ProblemData &problem_data, casadi::Function &upper_bound, casadi::Function &lower_bound);
            
        double mu_;
        ApproximationOrder approximation_order_;

    };            
}
}
}