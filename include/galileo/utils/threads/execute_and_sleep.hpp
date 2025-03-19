#ifndef __galileo_utils_threads_thread_pool_hpp__
#define __galileo_utils_threads_thread_pool_hpp__

#include <chrono>
#include <thread>

namespace galileo
{

    /**
     * Helper function to execute a given function and sleep for the remainder of the specified timing.
     * Will not interrupt the function if it is too slow.
     *
     * @tparam Functor : function type
     * @param f : callable object to execute
     * @param frequency : the frequency the function should run at.
     */
    template <typename Functor>
    void executeAndSleep(Functor f, double frequency)
    {
        using clock = std::chrono::high_resolution_clock;
        const auto start = clock::now();

        // Execute wrapped function
        f();

        // Compute desired duration rounded to clock decimation
        const std::chrono::duration<double> desiredDuration(1.0 / frequency);
        const auto dt = std::chrono::duration_cast<clock::duration>(desiredDuration);

        // Sleep
        const auto sleepTill = start + dt;
        std::this_thread::sleep_until(sleepTill);
    }

    /**
     * Helper function to execute a given function at a given rate while a condition is true.
     * Will not interrupt the function if it is too slow.
     *
     * @tparam Functor : function type
     * @param f : callable object to execute.
     * @param condition : condition checked every loop, will exit when this function returns false.
     * @param frequency : the frequency the function should run at.
     */
    template <typename Functor1, typename Functor2>
    void executeAtRate(Functor1 f, Functor2 condition, double frequency)
    {
        using clock = std::chrono::high_resolution_clock;

        // Compute desired duration rounded to clock decimation
        const std::chrono::duration<double> desiredDuration(1.0 / frequency);
        const auto dt = std::chrono::duration_cast<clock::duration>(desiredDuration);

        // Initialize timing
        const auto start = clock::now();
        auto sleepTill = start + dt;

        // Execution loop
        while (condition())
        {
            // Execute wrapped function
            f();

            // Sleep
            std::this_thread::sleep_until(sleepTill);
            sleepTill += dt;
        }
    }

} // namespace galileo

#endif // __galileo_utils_threads_thread_pool_hpp__
