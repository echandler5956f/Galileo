#ifndef __galileo_macro_hpp_
#define __galileo_macro_hpp_

#define GALILEO_WORLD_VERSION 2
#define GALILEO_MAJOR_VERSION 0
#define GALILEO_MINOR_VERSION 00

#define GALILEO_VERSION_AT_LEAST(x, y, z) (GALILEO_WORLD_VERSION > x || (GALILEO_WORLD_VERSION >= x &&                                \
                                                                         (GALILEO_MAJOR_VERSION > y || (GALILEO_MAJOR_VERSION >= y && \
                                                                                                        GALILEO_MINOR_VERSION >= z))))

#define FORWARD_ACCESSOR(ReturnType, accessor_name)  \
    /* lvalue-qualified overload */                  \
    ReturnType &accessor_name()                      \
    {                                                \
        return derived().accessor_name##_accessor(); \
    }                                                \
    /* const-lvalue-qualified overload */            \
    const ReturnType &accessor_name() const          \
    {                                                \
        return derived().accessor_name##_accessor(); \
    }

#define DEFAULT_ACCESSOR(ReturnType, accessor_name)    \
    ReturnType &accessor_name##_accessor()             \
    {                                                  \
        return accessor_name;                          \
    }                                                  \
    const ReturnType &accessor_name##_accessor() const \
    {                                                  \
        return accessor_name;                          \
    }

#define GENERIC_ACCESSOR(ReturnType, accessor_name)    \
    ReturnType &accessor_name##_accessor()             \
    {                                                  \
        return accessor_name();                        \
    }                                                  \
    const ReturnType &accessor_name##_accessor() const \
    {                                                  \
        return accessor_name();                        \
    }

#define GALILEO_DEFAULT_CONSTRUCTOR(X) \
    X() = default;                     \
    ~X() = default;                    \
    X(const X &) = default;            \
    X(X &&) = default;

#endif // __galileo_macro_hpp_