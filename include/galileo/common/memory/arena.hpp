#ifndef __galileo_common_memory_arena_hpp__
#define __galileo_common_memory_arena_hpp__

#include "galileo/fwd.hpp"

#include <cstddef>
#include <cstdint>
#include <memory_resource>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace galileo
{

    namespace detail
    {
        inline constexpr std::size_t round_up(std::size_t n, std::size_t multiple) noexcept
        {
            return (n + multiple - 1u) / multiple * multiple;
        }

        inline std::byte *align_ptr(std::byte *p, std::size_t align) noexcept
        {
            // align is a power of two in PMR.
            auto addr = reinterpret_cast<std::uintptr_t>(p);
            auto mis = addr & (align - 1u);
            if (mis == 0) return p;
            return p + (align - mis);
        }

    } // namespace detail

    // -------------------------------------------------------------------------------------------------
    // MemoryArena: a strictly-aligned monotonic arena memory_resource
    // -------------------------------------------------------------------------------------------------
    class MemoryArena : public std::pmr::memory_resource
    {
    public:
        // Construct with an external buffer. The arena will start by carving allocations
        // out of that buffer and then fall back to upstream for additional, aligned blocks.
        MemoryArena(void *buffer,
                    std::size_t bytes,
                    std::pmr::memory_resource *upstream = std::pmr::get_default_resource(),
                    std::size_t block_growth = 64 * 1024,
                    std::size_t min_alignment = detail::k_default_align_bytes) noexcept
            : upstream_(upstream),
              min_align_(min_alignment < alignof(std::max_align_t) ? alignof(std::max_align_t) : min_alignment),
              growth_(block_growth)
        {
            if (buffer && bytes)
            {
                blocks_.push_back(Block{static_cast<std::byte *>(buffer),
                                        bytes,
                                        0u,
                                        /*own*/ false,
                                        /*align*/ min_align_});
            }
        }

        // Construct with no external buffer; first block will be allocated on demand from upstream.
        explicit MemoryArena(std::pmr::memory_resource *upstream,
                             std::size_t initial_block = 64 * 1024,
                             std::size_t min_alignment = detail::k_default_align_bytes,
                             std::size_t block_growth = 64 * 1024) noexcept
            : upstream_(upstream),
              min_align_(min_alignment < alignof(std::max_align_t) ? alignof(std::max_align_t) : min_alignment),
              growth_(block_growth),
              default_block_(initial_block)
        {
        }

        ~MemoryArena() override { release(); }

        // Drop all allocated memory at once.
        void release() noexcept
        {
            for (auto &b : blocks_)
            {
                if (b.own && b.begin)
                {
                    upstream_->deallocate(b.begin, b.capacity, b.align ? b.align : min_align_);
                }
            }
            blocks_.clear();
        }

    protected:
        void *do_allocate(std::size_t bytes, std::size_t alignment) override
        {
            const std::size_t A = alignment > min_align_ ? alignment : min_align_;

            // Fast path: try last block
            if (!blocks_.empty())
            {
                auto &b = blocks_.back();
                std::byte *base = detail::align_ptr(b.begin + b.offset, A);
                if (base + bytes <= b.begin + b.capacity)
                {
                    b.offset = static_cast<std::size_t>(base - b.begin) + bytes;
                    return base;
                }
            }

            // Need a new block
            const std::size_t want = detail::round_up(bytes + A, A);
            std::size_t block_size = blocks_.empty() ? (default_block_ ? std::max(default_block_, want) : want)
                                                     : std::max(blocks_.back().capacity * 2u, std::max(growth_, want));

            void *raw = upstream_->allocate(block_size, A);
            blocks_.push_back(Block{static_cast<std::byte *>(raw),
                                    block_size,
                                    0u,
                                    /*own*/ true,
                                    /*align*/ A});

            auto &b = blocks_.back();
            std::byte *base = detail::align_ptr(b.begin, A);
            b.offset = static_cast<std::size_t>(base - b.begin) + bytes;
            return base;
        }

        void do_deallocate(void * /*p*/, std::size_t /*bytes*/, std::size_t /*alignment*/) override
        {
            // Monotonic resource: individual deallocations are no-ops.
        }

        bool do_is_equal(const std::pmr::memory_resource &other) const noexcept override { return this == &other; }

    private:
        struct Block
        {
            std::byte *begin{};
            std::size_t capacity{};
            std::size_t offset{};
            bool own{};
            std::size_t align{}; // alignment used to allocate this block from upstream
        }; // struct Block

        std::pmr::memory_resource *upstream_ = std::pmr::get_default_resource();
        std::size_t min_align_{detail::k_default_align_bytes};
        std::size_t growth_{64 * 1024};
        std::size_t default_block_{64 * 1024};
        std::vector<Block> blocks_{};
    }; // class MemoryArena

} // namespace galileo

#endif // __galileo_common_memory_arena_hpp__
