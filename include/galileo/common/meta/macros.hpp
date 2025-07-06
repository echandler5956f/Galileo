#ifndef __galileo_common_meta_macros_hpp__
#define __galileo_common_meta_macros_hpp__

#define GALILEO_WORLD_VERSION 2
#define GALILEO_MAJOR_VERSION 0
#define GALILEO_MINOR_VERSION 00

#define GALILEO_VERSION_AT_LEAST(x, y, z) (GALILEO_WORLD_VERSION > x || (GALILEO_WORLD_VERSION >= x &&                                \
                                                                         (GALILEO_MAJOR_VERSION > y || (GALILEO_MAJOR_VERSION >= y && \
                                                                                                        GALILEO_MINOR_VERSION >= z))))

#define FORWARD_ACCESSOR(ReturnType, accessor_name)        \
    /* lvalue-qualified overload */                        \
    ReturnType &accessor_name()                            \
    {                                                      \
        return this->derived().accessor_name##_accessor(); \
    }                                                      \
    /* const-lvalue-qualified overload */                  \
    const ReturnType &accessor_name() const                \
    {                                                      \
        return this->derived().accessor_name##_accessor(); \
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

#endif // __galileo_common_meta_macros_hpp__
