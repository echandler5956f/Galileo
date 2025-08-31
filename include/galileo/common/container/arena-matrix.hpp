#ifndef __galileo_common_container_arena_matrix_hpp__
#define __galileo_common_container_arena_matrix_hpp__

#include <memory>
#include <type_traits>
#include <utility>

#include <Eigen/Core>

#include "galileo/common/meta/concepts.hpp"
#include "galileo/common/meta/dimension.hpp"
#include "galileo/common/meta/eigen.hpp"
#include "galileo/common/memory/arena.hpp"

namespace galileo
{

    // template <typename PlainLike>
    //     requires IsEigenDenseBase<PlainLike>
    // class ArenaMatrixTpl : public Eigen::Map<std::decay_t<PlainLike>, Eigen::Aligned>
    // {
    //     using Plain = std::decay_t<PlainLike>;
    //     using Base = Eigen::Map<Plain, Eigen::Aligned>;

    // public:
    //     using Scalar = typename Plain::Scalar;
    //     using Index = Eigen::Index;

    //     using DimRows = DimensionTpl<Plain::RowsAtCompileTime>;
    //     using DimCols = DimensionTpl<Plain::ColsAtCompileTime>;
    //     using DimMaxRows = DimensionTpl<Plain::MaxRowsAtCompileTime>;
    //     using DimMaxCols = DimensionTpl<Plain::MaxColsAtCompileTime>;

    //     // Lifetime management for non-trivial scalars
    //     static constexpr bool k_trivially_default_constructible = std::is_trivially_default_constructible_v<Scalar>;
    //     static constexpr bool k_trivially_destructible = std::is_trivially_destructible_v<Scalar>;
    //     static constexpr bool k_manage_lifetime = !(k_trivially_default_constructible && k_trivially_destructible);

    //     static constexpr std::size_t required_alignment() noexcept
    //     {
    //         constexpr std::size_t s_align = alignof(Scalar);
    //         return s_align > detail::k_default_align_bytes ? s_align : detail::k_default_align_bytes;
    //     }

    //     ArenaMatrixTpl(MemoryArena &arena,
    //                    DimRows rows = DimRows{},
    //                    DimCols cols = DimCols{},
    //                    DimMaxRows cap_rows = DimMaxRows{},
    //                    DimMaxCols cap_cols = DimMaxCols{})
    //         : Base(nullptr, rows.value(), cols.value())
    //     {
    //         construct_storage(arena, rows, cols, cap_rows, cap_cols);
    //     }

    //     ArenaMatrixTpl(MemoryArena &arena,
    //                    int rows = DimRows::Value,
    //                    int cols = DimCols::Value,
    //                    int cap_rows = DimMaxRows::Value,
    //                    int cap_cols = DimMaxCols::Value)
    //         : Base(nullptr, DimRows(rows).value(), DimCols(cols).value())
    //     {
    //         construct_storage(arena, DimRows(rows), DimCols(cols), DimMaxRows(cap_rows), DimMaxCols(cap_cols));
    //     }

    //     // Non-copyable (owning view)
    //     ArenaMatrixTpl(const ArenaMatrixTpl &) = delete;
    //     ArenaMatrixTpl &operator=(const ArenaMatrixTpl &) = delete;

    //     // Movable
    //     ArenaMatrixTpl(ArenaMatrixTpl &&other) noexcept
    //         : Base(other.data(), other.rows(), other.cols()),
    //           arena_(other.arena_),
    //           ptr_(other.ptr_),
    //           cap_rows_(other.cap_rows_),
    //           cap_cols_(other.cap_cols_)
    //     {
    //         other.remap_null();
    //     }

    //     ArenaMatrixTpl &operator=(ArenaMatrixTpl &&other) noexcept
    //     {
    //         if (this != &other)
    //         {
    //             // Destroy our current elements if needed
    //             if constexpr (k_manage_lifetime)
    //             {
    //                 if (ptr_)
    //                     std::destroy_n(ptr_, static_cast<std::size_t>(cap_rows_) * static_cast<std::size_t>(cap_cols_));
    //             }
    //             // Map to other's region
    //             new (this) Base(other.data(), other.rows(), other.cols());
    //             arena_ = other.arena_;
    //             ptr_ = other.ptr_;
    //             cap_rows_ = other.cap_rows_;
    //             cap_cols_ = other.cap_cols_;
    //             other.remap_null();
    //         }
    //         return *this;
    //     }

    //     ~ArenaMatrixTpl()
    //     {
    //         if constexpr (k_manage_lifetime)
    //         {
    //             if (ptr_)
    //             {
    //                 std::destroy_n(ptr_, static_cast<std::size_t>(cap_rows_) * static_cast<std::size_t>(cap_cols_));
    //             }
    //         }
    //     }

    //     // Capacity/introspection
    //     [[nodiscard]] Index capacity_rows() const noexcept { return cap_rows_; }
    //     [[nodiscard]] Index capacity_cols() const noexcept { return cap_cols_; }
    //     [[nodiscard]] std::size_t capacity_bytes() const noexcept
    //     {
    //         return static_cast<std::size_t>(cap_rows_) * static_cast<std::size_t>(cap_cols_) * sizeof(Scalar);
    //     }

    //     // In-capacity resize by remapping the base Map.
    //     void resize(DimRows new_rows = DimRows{}, DimCols new_cols = DimCols{})
    //     {
    //         if constexpr (DimRows::IsDynamic)
    //         {
    //             GALILEO_ASSERT(new_rows <= cap_rows_, "Rows exceed capacity");
    //         }
    //         if constexpr (DimCols::IsDynamic)
    //         {
    //             GALILEO_ASSERT(new_cols <= cap_cols_, "Cols exceed capacity");
    //         }
    //         new (this) Base(ptr_, new_rows.value(), new_cols.value());
    //     }

    // protected:
    //     void construct_storage(MemoryArena &arena,
    //                            DimRows rows = DimRows{},
    //                            DimCols cols = DimCols{},
    //                            DimMaxRows cap_rows = DimMaxRows{},
    //                            DimMaxCols cap_cols = DimMaxCols{})
    //     {
    //         // Validate fixed dims and clamp runtime ones
    //         if constexpr (DimRows::IsDynamic)
    //         {
    //             GALILEO_ASSERT(rows <= cap_rows, "Rows exceed capacity");
    //         }
    //         if constexpr (DimCols::IsDynamic)
    //         {
    //             GALILEO_ASSERT(cols <= cap_cols, "Cols exceed capacity");
    //         }

    //         const std::size_t total_elems = static_cast<std::size_t>(cap_rows) * static_cast<std::size_t>(cap_cols);
    //         const std::size_t total_bytes = total_elems * sizeof(Scalar);

    //         void *raw = arena.allocate(total_bytes, required_alignment());
    //         ptr_ = static_cast<Scalar *>(raw);

    //         if constexpr (k_manage_lifetime)
    //         {
    //             std::uninitialized_default_construct_n(ptr_, total_elems);
    //         }

    //         // Initialize the base Map
    //         new (static_cast<Base *>(this)) Base(ptr_, rows, cols);

    //         arena_ = &arena;
    //         cap_rows_ = cap_rows.value();
    //         cap_cols_ = cap_cols.value();
    //     }

    //     void remap_null() noexcept
    //     {
    //         // Rebind to a harmless null/zero-sized map and clear ownership flags
    //         new (static_cast<Base *>(this)) Base(static_cast<Scalar *>(nullptr),
    //                                              DimRows::IsDynamic ? Index(0) : Index(DimRows::Value),
    //                                              DimCols::IsDynamic ? Index(0) : Index(DimCols::Value));
    //         arena_ = nullptr;
    //         ptr_ = nullptr;
    //         cap_rows_ = 0;
    //         cap_cols_ = 0;
    //     }

    //     MemoryArena *arena_{};
    //     Scalar *ptr_{};
    //     Index cap_rows_{};
    //     Index cap_cols_{};
    // };

    // -------------------------------------------------------------------------------------------------
    // ArenaMatrixTpl: Eigen::Map backed by a std::pmr::memory_resource (Stan-style semantics).
    //
    //  - Default ctor maps to nullptr with 0 (or fixed) dims.
    //  - (arena, rows, cols) allocates rows*cols scalars in the arena and maps to them.
    //  - (arena, size) vector-only convenience ctor.
    //  - Expression ctor: allocates and assigns.
    //  - operator=(Eigen dense): if shape differs or null -> allocate & remap, then assign; else in-place assign.
    //  - Copy ctor / copy assign: shallow aliasing (no allocation), like Stan.
    //
    //  Notes:
    //    * Memory is never deallocated individually (monotonic arena).
    //    * For non-trivial Scalar types, this mirrors Stan’s semantics (no per-element lifetime mgmt).
    // -------------------------------------------------------------------------------------------------
    template <typename PlainLike>
        requires IsEigenDenseBase<PlainLike>
    class ArenaMatrixTpl : public Eigen::Map<std::decay_t<PlainLike>>
    {
        using Plain = std::decay_t<PlainLike>;
        using Base = Eigen::Map<Plain>;

    public:
        using Scalar = typename Plain::Scalar;
        using Index = Eigen::Index;

        static constexpr int RowsAtCompileTime = Plain::RowsAtCompileTime;
        static constexpr int ColsAtCompileTime = Plain::ColsAtCompileTime;

        static constexpr std::size_t required_alignment() noexcept
        {
            constexpr std::size_t s_align = alignof(Scalar);
            return s_align > detail::k_default_align_bytes ? s_align : detail::k_default_align_bytes;
        }

        // Default constructor: null map with 0 (or fixed) dims; uses default PMR for future allocations.
        ArenaMatrixTpl()
            : Base(static_cast<Scalar *>(nullptr),
                   RowsAtCompileTime == Eigen::Dynamic ? Index(0) : Index(RowsAtCompileTime),
                   ColsAtCompileTime == Eigen::Dynamic ? Index(0) : Index(ColsAtCompileTime)),
              mr_(std::pmr::get_default_resource())
        {
        }

        // Construct with explicit arena and shape (matrix/array).
        ArenaMatrixTpl(std::pmr::memory_resource &arena, Index rows, Index cols)
            : Base(static_cast<Scalar *>(
                       arena.allocate(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols) * sizeof(Scalar),
                                      required_alignment())),
                   rows,
                   cols),
              mr_(&arena)
        {
        }

        // Construct with explicit arena and size (vector-only).
        explicit ArenaMatrixTpl(std::pmr::memory_resource &arena, Index size)
            : Base(static_cast<Scalar *>(
                       arena.allocate(static_cast<std::size_t>(size) * sizeof(Scalar), required_alignment())),
                   size),
              mr_(&arena)
        {
        }

    private:
        template <typename T>
        static Index get_rows(const T &x)
        {
            return (RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1) ||
                    (ColsAtCompileTime == 1 && T::RowsAtCompileTime == 1)
                ? x.cols()
                : x.rows();
        }
        template <typename T>
        static Index get_cols(const T &x)
        {
            return (RowsAtCompileTime == 1 && T::ColsAtCompileTime == 1) ||
                    (ColsAtCompileTime == 1 && T::RowsAtCompileTime == 1)
                ? x.rows()
                : x.cols();
        }

    public:
        // Expression constructor with explicit arena.
        template <typename T>
            requires IsEigenDenseBase<T>
        ArenaMatrixTpl(std::pmr::memory_resource &arena, const T &other) // NOLINT
            : Base(static_cast<Scalar *>(
                       arena.allocate(static_cast<std::size_t>(other.size()) * sizeof(Scalar), required_alignment())),
                   get_rows(other),
                   get_cols(other)),
              mr_(&arena)
        {
            Base::operator=(other);
        }

        // Expression constructor using default PMR.
        template <typename T>
            requires IsEigenDenseBase<T>
        explicit ArenaMatrixTpl(const T &other) // NOLINT
            : ArenaMatrixTpl(*std::pmr::get_default_resource(), other)
        {
        }

        // Construct from an existing Map (assumed arena-backed); shallow alias.
        explicit ArenaMatrixTpl(const Base &other) // NOLINT
            : Base(other), mr_(std::pmr::get_default_resource())
        {
        }

        // Copy constructor: shallow alias (no allocation).
        ArenaMatrixTpl(const ArenaMatrixTpl &other)
            : Base(const_cast<Scalar *>(other.data()), other.rows(), other.cols()), mr_(other.mr_)
        {
        }

        // Bring Eigen::Map assignment operators (scalar assign, etc.) into scope.
        using Base::operator=;

        // Copy assignment: shallow alias (no allocation).
        ArenaMatrixTpl &operator=(const ArenaMatrixTpl &other)
        {
            new (static_cast<Base *>(this)) Base(const_cast<Scalar *>(other.data()), other.rows(), other.cols());
            mr_ = other.mr_;
            return *this;
        }

        // Assignment from any Eigen dense expression:
        //  - If current shape differs (or null), allocate new memory in our arena and remap, then assign.
        //  - If shape matches, perform in-place assignment (no allocation).
        template <typename T>
            requires IsEigenDenseBase<T>
        ArenaMatrixTpl &operator=(const T &other)
        {
            const Index r = get_rows(other);
            const Index c = get_cols(other);

            if (this->data() == nullptr || this->rows() != r || this->cols() != c)
            {
                Scalar *p = static_cast<Scalar *>(
                    mr_->allocate(static_cast<std::size_t>(other.size()) * sizeof(Scalar), required_alignment()));
                new (static_cast<Base *>(this)) Base(p, r, c);
            }
            Base::operator=(other);
            return *this;
        }

        // Force a hard copy into current storage (no rebind, shapes must match).
        template <typename T>
        void deep_copy(const T &x)
        {
            Base::operator=(x);
        }

    private:
        std::pmr::memory_resource *mr_{};
    };

} // namespace galileo

#endif // __galileo_common_container_arena_matrix_hpp__
