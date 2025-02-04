#ifndef __galileo_simulator_simulator_base_hpp__
#define __galileo_simulator_simulator_base_hpp__

#include "galileo/simulator/fwd.hpp"

#include <functional>
#include <bits/std_function.h>

namespace galileo
{
    namespace simulator
    {

        class SimulatorBase
        {
        public:
            SimulatorBase(std::function<void(void)> loop_func, std::function<void(void)> exit_func = []() {}) : loop_func_(loop_func), exit_func_(exit_func) {}
            SimulatorBase(const SimulatorBase &) = delete;
            SimulatorBase &operator=(const SimulatorBase &) = delete;
            virtual ~SimulatorBase() = default;

            virtual void Initialize(const std::string &model_path)
            {
            }

            virtual void Configure(void *config)
            {
            }

            virtual void Finalize()
            {
            }

            void Loop()
            {
                while (true)
                {
                    LoopPrologue();
                    loop_func_();
                    LoopEpilogue();
                }
            }

            virtual void Exit()
            {
                exit_func_();
            }

        protected:
            virtual void LoopPrologue()
            {
            }

            virtual void LoopEpilogue()
            {
            }

            std::function<void(void)> loop_func_;
            std::function<void(void)> exit_func_;

        }; // class SimulatorBase

    } // namespace simulator

} // namespace galileo

#endif // __galileo_simulator_simulator_base_hpp__