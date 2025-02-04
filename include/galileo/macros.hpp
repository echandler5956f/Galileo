#ifndef __galileo_macro_hpp_
#define __galileo_macro_hpp_

#define GALILEO_WORLD_VERSION 2
#define GALILEO_MAJOR_VERSION 0
#define GALILEO_MINOR_VERSION 00

#define GALILEO_VERSION_AT_LEAST(x, y, z) (GALILEO_WORLD_VERSION > x || (GALILEO_WORLD_VERSION >= x &&                                \
                                                                         (GALILEO_MAJOR_VERSION > y || (GALILEO_MAJOR_VERSION >= y && \
                                                                                                        GALILEO_MINOR_VERSION >= z))))

#define FORWARD_GETTER(getter_name)       \
    /* lvalue-qualified overload */       \
    decltype(auto) getter_name() &        \
    {                                     \
        return derived().getter_name();   \
    }                                     \
    /* const-lvalue-qualified overload */ \
    decltype(auto) getter_name() const &  \
    {                                     \
        return derived().getter_name();   \
    }

#define GALILEO_DEFAULT_CONSTRUCTOR(X) \
    X() = default;                     \
    ~X() = default;                    \
    X(const X &) = default;            \
    X(X &&) = default;

#endif // __galileo_macro_hpp_