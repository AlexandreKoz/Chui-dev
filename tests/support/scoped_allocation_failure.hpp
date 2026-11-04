#pragma once

// Include in exactly one translation unit of a standalone test executable.
// Failure injection is thread-local, one-shot, and disabled while constructing
// exception diagnostics. No production hook or environment variable is added.
#include <cstddef>
#include <cstdlib>
#include <new>

namespace cosmosim::test_support {
inline thread_local std::ptrdiff_t allocations_until_failure = -1;
inline thread_local std::size_t failed_allocation_bytes = 0U;
inline thread_local std::size_t selected_allocation_bytes = 0U;

inline void allocationPoint(std::size_t bytes) {
  if (allocations_until_failure < 0) return;
  if (selected_allocation_bytes != 0U && selected_allocation_bytes != bytes) return;
  if (allocations_until_failure-- == 0) {
    allocations_until_failure = -1;
    failed_allocation_bytes = bytes;
    throw std::bad_alloc();
  }
}

class ScopedAllocationFailure {
 public:
  explicit ScopedAllocationFailure(std::ptrdiff_t count, std::size_t bytes = 0U) noexcept {
    failed_allocation_bytes = 0U;
    selected_allocation_bytes = bytes;
    allocations_until_failure = count;
  }
  ~ScopedAllocationFailure() { allocations_until_failure = -1; }
  ScopedAllocationFailure(const ScopedAllocationFailure&) = delete;
  ScopedAllocationFailure& operator=(const ScopedAllocationFailure&) = delete;
};
}  // namespace cosmosim::test_support

void* operator new(std::size_t bytes) {
  cosmosim::test_support::allocationPoint(bytes);
  if (void* p = std::malloc(bytes == 0U ? 1U : bytes)) return p;
  throw std::bad_alloc();
}
void* operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

void* operator new(std::size_t bytes, std::align_val_t alignment) {
  cosmosim::test_support::allocationPoint(bytes);
  const auto a = static_cast<std::size_t>(alignment);
  // aligned_alloc requires a multiple of alignment. Avoid size overflow.
  const auto padding = (a - bytes % a) % a;
  if (bytes > static_cast<std::size_t>(-1) - padding) throw std::bad_alloc();
  const auto size = bytes == 0U ? a : bytes + padding;
  if (void* p = std::aligned_alloc(a, size)) return p;
  throw std::bad_alloc();
}
void* operator new[](std::size_t bytes, std::align_val_t a) { return ::operator new(bytes, a); }
void operator delete(void* p, std::align_val_t) noexcept { std::free(p); }
void operator delete[](void* p, std::align_val_t) noexcept { std::free(p); }
void operator delete(void* p, std::size_t, std::align_val_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t, std::align_val_t) noexcept { std::free(p); }
