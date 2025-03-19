#ifndef __galileo_utils_threads_set_thread_priority_hpp__
#define __galileo_utils_threads_set_thread_priority_hpp__

#include <pthread.h>
#include <iostream>
#include <thread>

#include <bits/types/struct_sched_param.h>

namespace galileo
{

    /**
     * Sets the priority of the input thread.
     *
     * @param priority: The priority of the thread from 0 (lowest) to 99 (highest)
     * @param thread: A reference to the tread.
     */
    inline void setThreadPriority(int priority, pthread_t thread)
    {
        sched_param sched{};
        sched.sched_priority = priority;

        if (priority != 0)
        {
            if (pthread_setschedparam(thread, SCHED_FIFO, &sched) != 0)
            {
                std::cerr << "WARNING: Failed to set threads priority (one possible reason could be "
                             "that the user and the group permissions are not set properly.)"
                          << std::endl;
            }
        }
    }

    /**
     * Sets the priority of the input thread.
     *
     * @param priority: The priority of the thread from 0 (lowest) to 99 (highest)
     * @param thread: A reference to the tread.
     */
    inline void setThreadPriority(int priority, std::thread &thread)
    {
        setThreadPriority(priority, thread.native_handle());
    }

    /**
     * Sets the priority of the thread this function is called from.
     *
     * @param priority: The priority of the thread from 0 (lowest) to 99 (highest)
     */
    inline void setThisThreadPriority(int priority)
    {
        setThreadPriority(priority, pthread_self());
    }

} // namespace galileo

#endif // __galileo_utils_threads_set_thread_priority_hpp__
