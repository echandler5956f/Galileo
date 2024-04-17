#include "galileo/opt/ProblemData.h"
#include "galileo/opt/TrajectoryOpt.h"
#include "galileo/opt/States.h"
#include "galileo/opt/PhaseSequence.h"
#include "galileo/tools/MeshcatInterface.h"
#include "galileo/tools/GNUPlotInterface.h"

using namespace galileo;
using namespace opt;
using namespace tools;

struct SimpleProblemData
{
    std::shared_ptr<BasicSequence> phase_sequence_;
    std::shared_ptr<BasicStates> states;
    casadi::SX x;
    casadi::SX u;
    casadi::SX t;
};

template <class ProblemData>
class SimpleConstraintBuilder : public ConstraintBuilder<ProblemData>
{
public:
    SimpleConstraintBuilder() : ConstraintBuilder<ProblemData>() {}

    void buildConstraint(const ProblemData &problem_data, int phase_index, ConstraintData &constraint_data) override
    {
    } 

    void createBounds(const ProblemData &problem_data, int phase_index, casadi::Function &lower_bound, casadi::Function &upper_bound) const
    {
    }

    void createFunction(const ProblemData &problem_data, int phase_index, casadi::Function &G) const
    {
    }
};

template <class ProblemData>
class SimpleDecisionDataBuilder : public DecisionDataBuilder<ProblemData>
{
public:
    SimpleDecisionDataBuilder() : DecisionDataBuilder<ProblemData>() {}

    void buildDecisionData(const ProblemData &problem_data, int phase_index, DecisionData &decision_data) override
    {
        decision_data.lower_bound = casadi::Function("lower_bound", {problem_data.simple_constraint_problem_data.t}, {vertcat(casadi::SXVector{-0.25, -casadi::inf}), -1});
        decision_data.upper_bound = casadi::Function("upper_bound", {problem_data.simple_constraint_problem_data.t}, {vertcat(casadi::SXVector{casadi::inf, casadi::inf}), 1});
    }
};

class RosenbrockProblemData
{
public:
    RosenbrockProblemData(std::shared_ptr<GeneralProblemData> gp_data_,
                          std::shared_ptr<BasicSequence> basic_sequence_,
                          std::shared_ptr<BasicStates> states_,
                          casadi::SX x,
                          casadi::SX u,
                          casadi::SX t)
    {
        this->gp_data = gp_data_;
        this->phase_sequence = basic_sequence_;
        this->states = states_;

        this->simple_constraint_problem_data.phase_sequence_ = basic_sequence_;
        this->simple_constraint_problem_data.states = states_;
        this->simple_constraint_problem_data.x = x;
        this->simple_constraint_problem_data.u = u;
        this->simple_constraint_problem_data.t = t;
    }
    std::shared_ptr<PhaseSequence<BasicMode>> phase_sequence;
    std::shared_ptr<GeneralProblemData> gp_data;
    std::shared_ptr<BasicStates> states;
    SimpleProblemData simple_constraint_problem_data;
};