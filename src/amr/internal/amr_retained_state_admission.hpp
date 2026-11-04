#pragma once

#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cosmosim/core/memory_accounting.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::amr::internal {

// A narrow AMR persistent-state replacement transaction. The old owner stays
// in the process baseline while the complete replacement is physically live.
// The caller must reserve before constructing it and supply actual capacities.
// No allocation or scientific mutation is permitted in the commit callback.
class AmrRetainedStateAdmission final {
 public:
  AmrRetainedStateAdmission(const core::SimulationState& state,
                            core::MemoryGovernor* governor,
                            std::uint64_t replacement_bound,
                            std::string_view owner)
      : m_governor(governor) {
    if (governor == nullptr) return;
    m_state_before = core::memoryReportBaselineOwnedBytes(
        core::collectSimulationMemoryReport(state));
    const auto snapshot = governor->snapshot();
    m_process_before = snapshot.baseline_owned_bytes;
    if (m_process_before < m_state_before) {
      throw std::logic_error("AMR retained replacement has stale process baseline");
    }
    m_reservation = governor->reserve(core::MemoryClass::kCanonicalPersistent,
                                      replacement_bound, owner);
    m_reservation.commit();
  }

  AmrRetainedStateAdmission(const AmrRetainedStateAdmission&) = delete;
  AmrRetainedStateAdmission& operator=(const AmrRetainedStateAdmission&) = delete;

  // In-place retirement needs no candidate allocation. Validate every byte
  // calculation before the nonthrowing erase, then reconcile freed capacity.
  template <class TCommit>
  void commitShrink(std::uint64_t old_owned, std::uint64_t new_owned,
                    TCommit&& commit_state) {
    static_assert(std::is_nothrow_invocable_v<TCommit>);
    if (m_consumed || new_owned > old_owned || m_reservation.bytes() != 0U) {
      throw std::logic_error("invalid AMR retained in-place shrink");
    }
    if (m_governor != nullptr) {
      if (old_owned > m_state_before || old_owned > m_process_before) {
        throw std::logic_error("AMR retained shrink has stale baseline");
      }
      const auto baseline_after = m_process_before - (old_owned - new_owned);
      std::forward<TCommit>(commit_state)();
      m_consumed = true;
      m_reservation.reconcileBaselineOwnedAndRelease(baseline_after);
    } else {
      std::forward<TCommit>(commit_state)();
      m_consumed = true;
    }
  }

  template <class TCommit>
  void commit(std::uint64_t old_owned, std::uint64_t new_owned,
              TCommit&& commit_state) {
    static_assert(std::is_nothrow_invocable_v<TCommit>,
                  "AMR retained commit must be nonthrowing");
    if (m_consumed) throw std::logic_error("AMR retained admission reused");
    if (m_governor != nullptr) {
      if (new_owned > m_reservation.bytes() || old_owned > m_state_before ||
          old_owned > m_process_before) {
        throw std::logic_error("AMR retained replacement capacity exceeds admission");
      }
      const std::uint64_t baseline_after = core::checkedMemoryBytesAdd(
          m_process_before - old_owned, new_owned, "AMR retained replacement baseline");
      // All potentially throwing arithmetic is completed before canonical state changes.
      std::forward<TCommit>(commit_state)();
      m_consumed = true;
      m_reservation.reconcileBaselineOwnedAndRelease(baseline_after);
    } else {
      std::forward<TCommit>(commit_state)();
      m_consumed = true;
    }
  }

 private:
  core::MemoryGovernor* m_governor = nullptr;
  core::MemoryReservation m_reservation;
  std::uint64_t m_process_before = 0U;
  std::uint64_t m_state_before = 0U;
  bool m_consumed = false;
};

}  // namespace cosmosim::amr::internal
