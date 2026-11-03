#pragma once

#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cosmosim/core/memory_governor.hpp"

namespace cosmosim::core {

// Part 1 stage-transaction substrate: reserve -> prepare -> validate ->
// commit / reject for one conservative operation.
//
// Guarantees (and their limits):
//  * reserve(): all simultaneously resident transient bytes are admitted
//    through the MemoryGovernor before any work. A refusal is RECORDED (not
//    thrown) so every participant still reaches validate().
//  * prepare(callable): the callable may write ONLY to staging storage owned
//    by the caller (solver scratch, staged conserved copies, event buffers).
//    Exceptions are recorded and classified: StageNumericalRetryRequest ->
//    kNumericalRetry, StageModelDowngradeRequest -> kModelDowngradeRequired
//    (reserved for Part 2 thermochemistry, audit F2: an unfinished chemistry
//    interval must be rejected, never accepted as partial output),
//    MemoryAdmissionError -> kResourceRefused, anything else -> kInvalidState.
//  * validate(local_ok): exactly ONE collective vote per transaction. A
//    rank that already failed votes "reject"; the vote is reduced with the
//    supplied StageCollectiveDecision so all ranks see the same decision and
//    none blocks in a later collective. The vote must be issued outside
//    OpenMP parallel regions and outside any critical section; it is the
//    only collective the transaction performs.
//  * commit(callable): runs only after a unanimous accept. The callable moves
//    staged results into authoritative state and must not fail (it may only
//    write into pre-admitted storage). Pending reservations are committed
//    first so accounting matches residency. If the callable throws anyway the
//    transaction reports kCommitFailed and rethrows: that is an invariant
//    violation, not a recoverable rejection.
//  * Destruction releases every reservation exactly once (RAII), including on
//    exception paths.
//  * A rejected transaction leaves authoritative state untouched by
//    construction (nothing but commit() may write it). It does NOT roll back
//    work that a caller performed outside prepare()/commit(); cross-module
//    rollback (source + regrid + multirate) is Part 2/3 scope.
//  * Numerical retry (smaller dt / batch) is distinguishable from a physical
//    model downgrade; the caller decides policy. Part 1 callers treat any
//    rejection as fatal for the step after reporting it.

class StageNumericalRetryRequest : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

class StageModelDowngradeRequest : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

enum class StageTransactionPhase : std::uint8_t {
  kOpen = 0,
  kPrepared = 1,
  kValidated = 2,
  kCommitted = 3,
  kRejected = 4,
  kCommitFailed = 5,
};

enum class StageRejectionKind : std::uint8_t {
  kNone = 0,
  kResourceRefused = 1,
  kNumericalRetry = 2,
  kModelDowngradeRequired = 3,
  kInvalidState = 4,
  // This rank accepted, but at least one peer rejected.
  kPeerRejected = 5,
};

[[nodiscard]] std::string_view stageRejectionKindLabel(StageRejectionKind kind) noexcept;
[[nodiscard]] std::string_view stageTransactionPhaseLabel(StageTransactionPhase phase) noexcept;

// Collective agreement hook. allAccept() must be called by every participant
// of the communicator at the same transaction boundary.
class StageCollectiveDecision {
 public:
  virtual ~StageCollectiveDecision() = default;
  [[nodiscard]] virtual bool allAccept(bool local_accept, std::string_view phase) const = 0;
};

// Serial / non-MPI decision: the local vote is the decision.
class LocalStageCollectiveDecision final : public StageCollectiveDecision {
 public:
  [[nodiscard]] bool allAccept(bool local_accept, std::string_view phase) const override {
    (void)phase;
    return local_accept;
  }
};

struct StageTransactionOutcome {
  StageTransactionPhase phase = StageTransactionPhase::kOpen;
  StageRejectionKind rejection = StageRejectionKind::kNone;
  std::string reason;
  std::uint64_t reserved_bytes = 0U;
  [[nodiscard]] bool committed() const noexcept { return phase == StageTransactionPhase::kCommitted; }
};

class StageTransaction final {
 public:
  // governor may be null (unlimited/unaccounted execution); collective may be
  // null (equivalent to LocalStageCollectiveDecision).
  StageTransaction(std::string_view owner, MemoryGovernor* governor,
                   const StageCollectiveDecision* collective = nullptr);
  ~StageTransaction() noexcept;

  StageTransaction(const StageTransaction&) = delete;
  StageTransaction& operator=(const StageTransaction&) = delete;
  StageTransaction(StageTransaction&&) = delete;
  StageTransaction& operator=(StageTransaction&&) = delete;

  // Admit transient bytes. Returns false (and records kResourceRefused) on a
  // governor refusal; never throws MemoryAdmissionError.
  bool reserve(MemoryClass memory_class, std::uint64_t bytes, std::string_view label);
  // Adopt an already admitted reservation (pending or committed).
  void adoptReservation(MemoryReservation reservation);

  template <class TCallable>
  bool prepare(TCallable&& callable) {
    requirePhase(StageTransactionPhase::kOpen, "prepare");
    if (m_outcome.rejection != StageRejectionKind::kNone) {
      return false;
    }
    try {
      std::forward<TCallable>(callable)();
      m_outcome.phase = StageTransactionPhase::kPrepared;
      return true;
    } catch (const StageNumericalRetryRequest& error) {
      recordRejection(StageRejectionKind::kNumericalRetry, error.what());
    } catch (const StageModelDowngradeRequest& error) {
      recordRejection(StageRejectionKind::kModelDowngradeRequired, error.what());
    } catch (const MemoryAdmissionError& error) {
      recordRejection(StageRejectionKind::kResourceRefused, error.what());
    } catch (const std::exception& error) {
      recordRejection(StageRejectionKind::kInvalidState, error.what());
    }
    return false;
  }

  // The single collective vote. local_ok=false records `kind`/`reason`.
  // Returns the communicator-wide decision.
  bool validate(bool local_ok, StageRejectionKind kind_if_rejected, std::string_view reason_if_rejected);

  template <class TCallable>
  StageTransactionOutcome commit(TCallable&& callable) {
    requirePhase(StageTransactionPhase::kValidated, "commit");
    commitPendingReservations();
    try {
      std::forward<TCallable>(callable)();
    } catch (...) {
      m_outcome.phase = StageTransactionPhase::kCommitFailed;
      m_outcome.reason = std::string(m_owner) + ": commit callable failed after collective acceptance";
      throw;
    }
    m_outcome.phase = StageTransactionPhase::kCommitted;
    return m_outcome;
  }

  // Explicit local rejection before validate() (still votes in validate()).
  void reject(StageRejectionKind kind, std::string_view reason);

  [[nodiscard]] const StageTransactionOutcome& outcome() const noexcept { return m_outcome; }
  [[nodiscard]] std::string_view owner() const noexcept { return m_owner; }

 private:
  void requirePhase(StageTransactionPhase expected, const char* operation) const;
  void recordRejection(StageRejectionKind kind, std::string_view reason);
  void commitPendingReservations();

  std::string m_owner;
  MemoryGovernor* m_governor = nullptr;
  const StageCollectiveDecision* m_collective = nullptr;
  std::vector<MemoryReservation> m_reservations;
  StageTransactionOutcome m_outcome;
  bool m_voted = false;
};

// Thrown by callers that convert a rejected transaction into a failed step.
class StageTransactionRejected : public std::runtime_error {
 public:
  StageTransactionRejected(StageRejectionKind kind, const std::string& message)
      : std::runtime_error(message), m_kind(kind) {}
  [[nodiscard]] StageRejectionKind kind() const noexcept { return m_kind; }

 private:
  StageRejectionKind m_kind = StageRejectionKind::kNone;
};

}  // namespace cosmosim::core
