#ifndef __galileo_utils_threads_thread_pool_hpp__
#define __galileo_utils_threads_thread_pool_hpp__

#include "galileo/utils/threads/synchronized.hpp"

namespace galileo
{
    /**
     * Wraps a value with a thread-safe buffer. A value remains constant and can be accessed through the get function until
     * updateFromBuffer is called.
     *
     * In the meantime, multiple threads can set new values to the buffer. The active value is not protected by a mutex, so
     * only one thread should access/modify the active value (i.e. not simultaneously calling get() and updateFromBuffer()).
     *
     * @tparam T : wrapped type
     */
    template <typename T>
    class BufferedValue
    {
    public:
        /**
         * Constructor initializes with a given value and an empty buffer.
         * @param value
         */
        explicit BufferedValue(T value) : activeValue_(std::move(value)), buffer_(nullptr) {};

        /** Read the currently active value. */
        const T &get() const { return activeValue_; }

        /** Read/write the currently active value. */
        T &get() { return activeValue_; }

        /** Copy a new value into the buffer. */
        void setBuffer(const T &value) { buffer_.reset(std::make_unique<T>(value)); }

        /** Move a new value into the buffer. */
        void setBuffer(T &&value) { buffer_.reset(std::make_unique<T>(std::move(value))); }

        /**
         * Replaces the active value with the value in the buffer.
         * The active value is not mutex protected so this method is NOT thread-safe w.r.t. get()
         * The buffer is mutex protected, so this method is thread-safe w.r.t. setBuffer()
         * @return True: the active value was updated, False: the active value was not updated.
         */
        bool updateFromBuffer()
        {
            // Read buffer with a pointer swap to minimize time under the lock. The swapped value will be null if there was no new value set.
            std::unique_ptr<T> updatedValuePtr(nullptr);
            buffer_.swap(updatedValuePtr);

            if (updatedValuePtr != nullptr)
            {
                activeValue_ = std::move(*updatedValuePtr);
                return true;
            }
            else
            {
                return false;
            }
        }

    private:
        T activeValue_;
        Synchronized<T> buffer_;

    }; // class BufferedValue

} // namespace galileo

#endif // __galileo_utils_threads_thread_pool_hpp__
