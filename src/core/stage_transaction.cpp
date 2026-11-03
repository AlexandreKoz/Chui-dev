#include "cosmosim/core/stage_transaction.hpp"

namespace cosmosim::core {

std::string_view stageRejectionKindLabel(StageRejectionKind kind) noexcept {
  switch (kind) {
    case StageRejectionKind::kNone:
      return "none";
    case StageRejectionKind::kResourceRefused:
      return "resource_refused";
    case StageRejectionKind::kNumericalRetry:
      return "numerical_retry";
    case StageRejectionKind::kModelDowngradeRequired:
      return "model_downgrade_required";
    case StageRejectionKind::kInvalidState:
      return "invalid_state";
    case StageRejectionKind::kPeerRejected:
      return "peer_rejected";
  }
  return "unknown";
}

std::string_view stageTransactionPhaseLabel(StageTransactionPhase phase) noexcept {
  switch (phase) {
    case StageTransactionPhase::kOpen:
      return "open";
    case StageTransactionPhase::kPrepared:
      return "prepared";
    case StageTransactionPhase::kValidated:
      return "validated";
    case StageTransactionPhase::kCommitted:
      return "committed";
    case StageTransactionPhase::kRejected:
      return "rejected";
    case StageTransactionPhase::kCommitFailed:
      return "commit_failed";
  }
  return "unknown";
}

StageTransaction::StageTransaction(
    std::string_view owner, MemoryGovernor* governor, const StageCollectiveDecision* collective)
    : m_owner(owner), m_governor(governor), m_collective(collective) {
  if (m_owner.empty()) {
    throw std::invalid_argument("StageTransaction requires a non-empty owner label");
  }
}

StageTransaction::~StageTransaction() noexcept {
  for (MemoryReservation& reservation : m_reservations) {
    reservation.release();
  }
}

bool StageTransaction::reserve(MemoryClass memory_class, std::uint64_t bytes, std::string_view label) {
  requirePhase(StageTransactionPhase::kOpen, "reserve");
  if (m_outcome.rejection != StageRejectionKind::kNone) {
    return false;
  }
  if (m_governor == nullptr || bytes == 0U) {
    return true;
  }
  try {
    m_reservations.push_back(m_governor->reserve(memory_class, bytes, label));
    m_outcome.reserved_bytes = checkedMemoryBytesAdd(
        m_outcome.reserved_bytes, bytes, "StageTransaction reserved bytes");
    return true;
  } catch (const MemoryAdmissionError& error) {
    recordRejection(StageRejectionKind::kResourceRefused, error.what());
  } catch (const std::overflow_error& error) {
    recordRejection(StageRejectionKind::kResourceRefused, error.what());
  }
  return false;
}

void StageTransaction::adoptReservation(MemoryReservation reservation) {
  if (!reservation.valid()) {
    return;
  }
  m_outcome.reserved_bytes = checkedMemoryBytesAdd(
      m_outcome.reserved_bytes, reservation.bytes(), "StageTransaction adopted bytes");
  m_reservations.push_back(std::move(reservation));
}

bool StageTransaction::validate(
    bool local_ok, StageRejectionKind kind_if_rejected, std::string_view reason_if_rejected) {
  if (m_voted) {
    throw std::logic_error(m_owner + ": StageTransaction::validate may vote only once");
  }
  m_voted = true;
  if (!local_ok && m_outcome.rejection == StageRejectionKind::kNone) {
    recordRejection(kind_if_rejected == StageRejectionKind::kNone
                        ? StageRejectionKind::kInvalidState : kind_if_rejected,
                    reason_if_rejected);
  }
  const bool local_accept = m_outcome.rejection == StageRejectionKind::kNone &&
      (m_outcome.phase == StageTransactionPhase::kPrepared || m_outcome.phase == StageTransactionPhase::kOpen);
  const bool all_accept = m_collective != nullptr
      ? m_collective->allAccept(local_accept, m_owner)
      : local_accept;
  if (!all_accept) {
    if (m_outcome.rejection == StageRejectionKind::kNone) {
      recordRejection(StageRejectionKind::kPeerRejected, "a peer rank rejected the transaction");
    }
    m_outcome.phase = StageTransactionPhase::kRejected;
    for (MemoryReservation& reservation : m_reservations) {
      reservation.release();
    }
    return false;
  }
  m_outcome.phase = StageTransactionPhase::kValidated;
  return true;
}

void StageTransaction::reject(StageRejectionKind kind, std::string_view reason) {
  if (m_voted) {
    throw std::logic_error(m_owner + ": StageTransaction::reject after the collective vote");
  }
  recordRejection(kind == StageRejectionKind::kNone ? StageRejectionKind::kInvalidState : kind, reason);
}

void StageTransaction::requirePhase(StageTransactionPhase expected, const char* operation) const {
  if (m_outcome.phase != expected) {
    throw std::logic_error(m_owner + ": StageTransaction::" + operation + " called in phase " +
                           std::string(stageTransactionPhaseLabel(m_outcome.phase)));
  }
}

void StageTransaction::recordRejection(StageRejectionKind kind, std::string_view reason) {
  if (m_outcome.rejection != StageRejectionKind::kNone) {
    return;  // keep the first cause
  }
  m_outcome.rejection = kind;
  m_outcome.reason = m_owner + ": " + std::string(reason);
}

void StageTransaction::commitPendingReservations() {
  for (MemoryReservation& reservation : m_reservations) {
    if (reservation.pending()) {
      reservation.commit();
    }
  }
}

}  // namespace cosmosim::core
