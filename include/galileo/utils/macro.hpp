#ifndef __galileo_utils_macro_hpp_
#define __galileo_utils_macro_hpp_

#include <stdexcept>

namespace galileo
{
    /**
     * @brief A namespaced runtime_error for easier catch
     */
    struct runtime_error : std::runtime_error
    {
        using std::runtime_error::runtime_error;
        using std::runtime_error::what;
    };

    /**
     * @brief A namespaced invalid_argument for easier catch
     */
    struct invalid_argument : std::invalid_argument
    {
        using std::invalid_argument::invalid_argument;
        using std::invalid_argument::what;
    };

    namespace detail
    {

        template <typename E, typename... Args>
        void
#if defined(__GNUC__) || defined(__clang__)
            __attribute__((noinline, cold, noreturn))
#elif defined(_MSC_VER)
            __declspec(noinline, noreturn)
#else
// nothing
#endif
            raise(Args &&...args)
        {
            throw E(std::forward<Args>(args)...);
        }

    } // namespace detail

} // namespace galileo

#ifdef NDEBUG
#ifndef GALILEO_NO_DEBUG
#define GALILEO_NO_DEBUG
#endif
#endif

// gcc expands __VA_ARGS___ before passing it into the macro.
// Visual Studio expands __VA_ARGS__ after passing it.
// This macro is a workaround to support both
#define __GALILEO_EXPAND(x) x

#define __GALILEO_THROW_EXCEPT(msg, except) \
    galileo::detail::raise<except>(msg);

#define __GALILEO_THROW(msg) \
    __GALILEO_THROW_EXCEPT(msg, galileo::runtime_error)

#define __GALILEO_GET_MACRO_2(_1, _2, NAME, ...) NAME

#define GALILEO_THROW(...)                            \
    __GALILEO_EXPAND(                                 \
        __GALILEO_GET_MACRO_2(__VA_ARGS__,            \
                              __GALILEO_THROW_EXCEPT, \
                              __GALILEO_THROW)(__VA_ARGS__))

#define __GALILEO_CHECK_MSG_EXCEPT(cond, msg, except) \
    if (!(cond))                                      \
    {                                                 \
        GALILEO_THROW(msg, except);                   \
    }

#define __GALILEO_CHECK_MSG(cond, msg) \
    __GALILEO_CHECK_MSG_EXCEPT(cond, msg, galileo::runtime_error)

#define __GALILEO_CHECK(cond)   \
    __GALILEO_CHECK_MSG_EXCEPT( \
        cond, "Condition: '" #cond "' failed!", galileo::runtime_error)

#define __GALILEO_GET_MACRO_3(_1, _2, _3, NAME, ...) NAME

#define GALILEO_CHECK(...)                                \
    __GALILEO_EXPAND(                                     \
        __GALILEO_GET_MACRO_3(__VA_ARGS__,                \
                              __GALILEO_CHECK_MSG_EXCEPT, \
                              __GALILEO_CHECK_MSG,        \
                              __GALILEO_CHECK)(__VA_ARGS__))

// Assertions cost run time and can be turned off.
// You can suppress GALILEO_ASSERT by defining
// GALILEO_NO_DEBUG before including manif headers.
// GALILEO_NO_DEBUG is undefined by default unless NDEBUG is defined.
#ifndef GALILEO_NO_DEBUG
#define GALILEO_ASSERT(...)                               \
    __GALILEO_EXPAND(                                     \
        __GALILEO_GET_MACRO_3(__VA_ARGS__,                \
                              __GALILEO_CHECK_MSG_EXCEPT, \
                              __GALILEO_CHECK_MSG,        \
                              __GALILEO_CHECK)(__VA_ARGS__))
#else
#define GALILEO_ASSERT(...) ((void)0)
#endif

#define FORWARD_GETTER(getter_name)                \
    /* lvalue-qualified overload */                \
    decltype(auto) getter_name() &                 \
    {                                              \
        return derived().getter_name();            \
    }                                              \
    /* const-lvalue-qualified overload */          \
    decltype(auto) getter_name() const &           \
    {                                              \
        return derived().getter_name();            \
    }                                              \

#define GALILEO_DEFAULT_CONSTRUCTOR(X) \
    X() = default;                     \
    ~X() = default;                    \
    X(const X &) = default;            \
    X(X &&) = default;

#endif // __galileo_utils_macro_hpp_