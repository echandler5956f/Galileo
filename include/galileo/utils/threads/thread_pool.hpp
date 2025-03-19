#ifndef __galileo_utils_threads_thread_pool_hpp__
#define __galileo_utils_threads_thread_pool_hpp__

#include <condition_variable>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace galileo
{

    /**
     * Thread pool class to execute tasks on multiple threads.
     */
    class ThreadPool
    {
    public:
        /**
         * Constructor
         *
         * @param [in] nThreads: Number of threads to launch in the pool
         * @param [in] priority: The worker thread priority
         */
        explicit ThreadPool(size_t nThreads = 1, int priority = 0);

        /**
         * Destructor
         */
        ~ThreadPool();

        /**
         * Run a task in another thread
         *
         * @tparam Functor: The task function
         * @param [in] taskFunction: The task function to run in the pool. It takes a thread worker index argument (between 0 and nThreads - 1),
                                     which can be used to index designated thread resources.
         * @return future object with taskFunction return value
         */
        template <typename Functor>
        std::future<typename std::result_of<Functor(int)>::type> run(Functor taskFunction);

        /**
         * Helper function to run a task N times parallel with the help of the pool.
         * - 1 task will run in the calling thread with ID = nThreads.
         * - N-1 tasks will run on the threadpool with ID in [0, nThreads-1].
         *
         * @note This is a blocking operation, returns when all tasks are completed.
         * @warning Calling runParallel(task, nThreads) does not guarantee that each task will be executed with a different workerIndex.
         *
         * @param [in] taskFunction: task function to run in the pool.
         * @param [in] N: number of times to run taskFunction in parallel.
         */
        void runParallel(std::function<void(int)> taskFunction, int N);

        /** Get the number of threads. */
        size_t numThreads() const { return workerThreads_.size(); }

    private:
        struct TaskBase; // forward declaration

        template <typename Functor>
        struct Task; // forward declaration

        /**
         * Thread worker loop
         *
         * @param [in] workerIndex: worker thread index
         */
        void worker(int workerIndex);

        /**
         * Run a task asynchronously in another thread
         *
         * @param [in] taskPtr: task object
         */
        void runTask(std::unique_ptr<TaskBase> taskPtr);

        bool stop_{false}; //!< flag telling all threads to stop, protected by taskQueueLock_

        std::queue<std::unique_ptr<TaskBase>> taskQueue_; // protected by taskQueueLock_
        std::condition_variable taskQueueCondition_;
        std::mutex taskQueueLock_;

        std::vector<std::thread> workerThreads_;

    }; // class ThreadPool

    /**
     * Task callback interface class.
     */
    struct ThreadPool::TaskBase
    {
        TaskBase() = default;
        virtual ~TaskBase() = default;
        virtual void operator()(int workerIndex) = 0;

    }; // struct ThreadPool::TaskBase

    /**
     * Task callback for a specific function return type.
     *
     * @tparam Functor: Type of callable task object (functor).
     */
    template <typename Functor>
    struct ThreadPool::Task final : public ThreadPool::TaskBase
    {
        explicit Task(Functor taskFunction) : packagedTask(std::move(taskFunction)) {}
        ~Task() override = default;
        void operator()(int workerIndex) override { packagedTask(workerIndex); }

        using ReturnType = typename std::result_of<Functor(int)>::type;
        std::packaged_task<ReturnType(int)> packagedTask;

    }; // struct ThreadPool::Task

    /**************************************************************************************************/
    /**************************************************************************************************/
    /**************************************************************************************************/
    template <typename Functor>
    std::future<typename std::result_of<Functor(int)>::type> ThreadPool::run(Functor taskFunction)
    {
        auto taskPtr = std::make_unique<Task<Functor>>(std::move(taskFunction));
        auto future = taskPtr->packagedTask.get_future();

        if (workerThreads_.empty())
        {
            // run on main thread
            taskPtr->operator()(0);
        }
        else
        {
            runTask(std::move(taskPtr));
        }

        return future;
    }

} // namespace galileo

#endif // __galileo_utils_threads_thread_pool_hpp__
