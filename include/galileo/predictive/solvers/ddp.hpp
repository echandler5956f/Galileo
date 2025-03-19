#ifndef __galileo_predictive_solvers_ddp_hpp__
#define __galileo_predictive_solvers_ddp_hpp__

#include "galileo/predictive/solvers/fwd.hpp"
#include "galileo/predictive/solvers/solver-base.hpp"

#include <vector>
#include <memory>
#include <limits>

namespace galileo
{

    namespace predictive
    {

        template <typename _VarScalar, typename _NumScalar, int _Options, template <typename, typename, int> class PhaseCollectionTpl>
        class SolverDDP; // forward declaration

        template <typename _VarScalar, typename _NumScalar, int _Options, template <typename, typename, int> class PhaseCollectionTpl>
        struct traits<SolverDDP<_VarScalar, _NumScalar, _Options, PhaseCollectionTpl>>
        {
            using SolverDerived = SolverDDP<_VarScalar, _NumScalar, _Options, PhaseCollectionTpl>;
            using VarScalar = _VarScalar;
            using NumScalar = _NumScalar;
            static constexpr int Options = _Options;

            using OptimalControlProblem_t = OptimalControlProblem<VarScalar, NumScalar, Options, PhaseCollectionTpl>;
            using VectorXv = Eigen::Matrix<VarScalar, Eigen::Dynamic, 1>;
            using VectorXn = Eigen::Matrix<NumScalar, Eigen::Dynamic, 1>;
        }; // struct traits

        template <typename _VarScalar, typename _NumScalar, int _Options, template <typename, typename, int> class PhaseCollectionTpl>
        class SolverDDP : SolverBase<SolverDDP<_VarScalar, _NumScalar, _Options, PhaseCollectionTpl>>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using SolverDerived = SolverDDP<_VarScalar, _NumScalar, _Options, PhaseCollectionTpl>;
            GALILEO_SOLVER_BASIC_TYPEDEF(SolverDerived);
            GALILEO_SOLVER_TYPEDEF(SolverDerived);

            using MatrixXv = Eigen::Matrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic>;
            using MatrixXn = Eigen::Matrix<NumScalar, Eigen::Dynamic, Eigen::Dynamic>;
            using MatrixXvRowMajor = Eigen::Matrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
            using MatrixXnsRowMajor = Eigen::Matrix<NumScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

            explicit SolverDDP(std::shared_ptr<ShootingProblem> problem);

            ~SolverDDP();

            bool solve(
                const std::vector<VectorXn> &init_xs,
                const std::vector<VectorXn> &init_us,
                const std::size_t maxiter = 100, const bool is_feasible = false,
                const NumScalar init_reg = NAN);

            void computeDirection(const bool recalc = true);

            NumScalar tryStep(const NumScalar steplength = 1);

            NumScalar stoppingCriteria();

            const Vector2ns &expectedImprovement();

            void resizeData();

            NumScalar calcDiff();

            void backwardPass();

            void forwardPass(const NumScalar stepLength);

            void computeActionValueFunction(
                const std::size_t t, const std::shared_ptr<ActionModelAbstract> &model,
                const std::shared_ptr<ActionDataAbstract> &data);

            void computeValueFunction(const std::size_t t, const std::shared_ptr<ActionModelAbstract> &model);

            void computeGains(const std::size_t t);

            void increaseRegularization();

            void decreaseRegularization();

            void allocateData();

            NumScalar get_reg_incfactor() const;

            NumScalar get_reg_decfactor() const;

            NumScalar get_reg_min() const;

            NumScalar get_reg_max() const;

            const std::vector<NumScalar> &get_alphas() const;

            NumScalar get_th_stepdec() const;

            NumScalar get_th_stepinc() const;

            NumScalar get_th_grad() const;

            const std::vector<MatrixXn> &get_Vxx() const;

            const std::vector<VectorXn> &get_Vx() const;

            const std::vector<MatrixXn> &get_Qxx() const;

            const std::vector<MatrixXn> &get_Qxu() const;

            const std::vector<MatrixXn> &get_Quu() const;

            const std::vector<VectorXn> &get_Qx() const;

            const std::vector<VectorXn> &get_Qu() const;

            const std::vector<MatrixXnsRowMajor> &get_K() const;

            const std::vector<VectorXn> &get_k() const;
            void set_reg_incfactor(const NumScalar reg_factor);

            void set_reg_decfactor(const NumScalar reg_factor);

            void set_reg_min(const NumScalar regmin);

            void set_reg_max(const NumScalar regmax);

            void set_alphas(const std::vector<NumScalar> &alphas);

            void set_th_stepdec(const NumScalar th_step);

            void set_th_stepinc(const NumScalar th_step);

            void set_th_grad(const NumScalar th_grad);

        protected:
            NumScalar reg_incfactor_;

            NumScalar reg_decfactor_;

            NumScalar reg_min_;
            NumScalar reg_max_;

            NumScalar cost_try_;
            std::vector<VectorXn> xs_try_;
            std::vector<VectorXn> us_try_;
            std::vector<VectorXn> dx_;

            std::vector<MatrixXn> Vxx_;
            MatrixXn Vxx_tmp_;
            std::vector<VectorXn> Vx_;
            std::vector<MatrixXn> Qxx_;
            std::vector<MatrixXn> Qxu_;
            std::vector<MatrixXn> Quu_;
            std::vector<VectorXn> Qx_;
            std::vector<VectorXn> Qu_;
            std::vector<MatrixXnsRowMajor> K_;
            std::vector<VectorXn> k_;

            VectorXn xnext_;
            MatrixXnsRowMajor FxTVxx_p_;

            std::vector<MatrixXnsRowMajor> FuTVxx_p_;

            VectorXn fTVxx_p_;

            std::vector<Eigen::LLT<MatrixXn>> Quu_llt_;
            std::vector<VectorXn> Quuk_;

            std::vector<NumScalar> alphas_;
            NumScalar th_grad_;

            NumScalar th_stepdec_;
            NumScalar th_stepinc_;

        }; // class SolverDDP

    } // namespace predictive

} // namespace galileo

// #include "galileo/predictive/solvers/ddp.hxx"

#endif // __galileo_predictive_solvers_ddp_hpp__