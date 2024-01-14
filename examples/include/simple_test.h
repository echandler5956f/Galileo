#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/opt/PseudospectralSegment.h"
#include <Eigen/Dense>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include "third-party/gnuplot-iostream/gnuplot-iostream.h"
#include "galileo/opt/BasicStates.h"

using namespace galileo;

struct RosenbrockProblemData
{
};

struct SimpleProblemData
{
    std::shared_ptr<opt::GeneralProblemData> gp_data;
    casadi::SX x;
    casadi::SX u;
    std::shared_ptr<opt::BasicStates> states;
    RosenbrockProblemData rosenbrock_problem_data;

    int num_knots;
};

template <class ProblemData>
class RosenbrockConstraintBuilder : public opt::ConstraintBuilder<ProblemData>
{
public:
    RosenbrockConstraintBuilder() : opt::ConstraintBuilder<ProblemData>() {}

    void CreateApplyAt(const ProblemData &problem_data, int knot_index, Eigen::VectorXi &apply_at) const
    {
        uint num_points = problem_data.num_knots;
        apply_at = Eigen::VectorXi::Constant(num_points, 1);
    }

    void CreateBounds(const ProblemData &problem_data, int knot_index, casadi::Function &upper_bound, casadi::Function &lower_bound) const
    {

    }

    void CreateFunction(const ProblemData &problem_data, int knot_index, casadi::Function &G) const
    {
        
    }
};