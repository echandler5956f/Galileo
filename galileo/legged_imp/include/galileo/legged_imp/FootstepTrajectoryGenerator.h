#include <eigen3/Eigen/Dense>
#include <iostream>

namespace galileo
{
    namespace model
    {
        namespace legged
        {
            class FootstepTrajectoryGenerator
            {
            public:
                FootstepTrajectoryGenerator() {}
                ~FootstepTrajectoryGenerator() = default;

                struct FootstepDefinition
                {
                    double h_start;
                    double h_end;
                    double h_max;
                };

                virtual Eigen::Vector2d calcFootstepZ(int collocation_point, int collocation_point_span, const FootstepDefinition &FS_def);
            };

            class HermiteFootstepGenerator : public FootstepTrajectoryGenerator
            {
                // Parameterizes footstep heights.
                // t = 0 is the start of the step
                // t = 1 is the end of the step
                // f(t) is the height of the footstep at time t

            public:
                HermiteFootstepGenerator() : FootstepTrajectoryGenerator() {}

                Eigen::Vector2d calcFootstepZ(int collocation_point, int collocation_point_span, const FootstepDefinition &FS_def)
                {
                    double t = std::static_cast<double>(collocation_point) /
                               std::static_cast<double>(collocation_point_span);

                    Eigen::Vector2d x0 = Eigen::Vector2d(FS_def.h_start, 0);
                    Eigen::Vector2d xf = Eigen::Vector2d(FS_def.h_end, 0);

                    if (t < 0.5)
                    {
                        xf(0) = FS_def.h_max;
                    }
                    else
                    {
                        x0(0) = FS_def.h_max;
                        t = t - 0.5;
                    }

                    return CubicHermiteSpline(t, x0, xf);
                }

            private:
                /**
                 * @brief Cubic Hermite Spline in one dimension
                 */
                Eigen::Vector2d CubicHermiteSpline(double t, Eigen::Vector2d x0, Eigen::Vector2d xf)
                {
                    Eigen::MatrixXd H0(2, 2);
                    H0(0, 0) = 2 * pow(t, 3) - 3 * pow(t, 2) + 1;
                    H0(0, 1) = pow(t, 3) - 2 * pow(t, 2) + t;
                    H0(1, 0) = 6 * pow(t, 2) - 6 * t;
                    H0(1, 1) = 3 * pow(t, 2) - 4 * t + 1;

                    Eigen::MatrixXd Hf(2, 2);
                    Hf(0, 0) = -2 * pow(t, 3) + 3 * pow(t, 2);
                    Hf(0, 1) = pow(t, 3) - pow(t, 2);
                    Hf(1, 0) = -6 * pow(t, 2) + 6 * t;
                    Hf(1, 1) = 3 * pow(t, 2) - 2 * t;

                    return H0 * x0 + Hf * xf;
                }
            };
        }
    }
}