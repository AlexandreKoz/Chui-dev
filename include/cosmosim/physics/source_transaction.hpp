#pragma once

// CHUI-RT-BH-AGN-REF-v1, shared contract B02: typed conservative source
// transactions for BH accretion/seeding/mergers, AGN injection and RT matter
// exchange (master design section 5; AGN companion section 3).
//
// Lifecycle (exactly once):
//   plan -> prepare (aggregate per recipient row, check identity, epoch,
//   finiteness, floors and positivity on a scratch image; NO state mutation)
//   -> validate (ledger closure of every plan against its declared source
//   increment) -> authorize (unanimous participant vote supplied by the
//   caller; a single-rank caller passes the local verdict) -> commit (apply
//   every aggregated row delta once, in ascending stable gas_cell_id order).
// A failed prepare/validate/authorize leaves canonical state untouched. A
// commit that is entered after authorization cannot fail for state reasons
// (prepare proved admissibility on the identical aggregated deltas); any
// exception thrown afterwards is a fatal internal error and the workflow must
// recover from the last committed checkpoint (no partial-continue path).
//
// Conservation convention (identical to MechanicalFeedbackDeposition::commit):
//   m' = m + dm,  p' = m v + dp,  E' = m (u + |v|^2/2) + dE_total,
//   v' = p'/m',   u' = (E' - |p'|^2/(2 m'))/m'.
// Withdrawal of a fraction f (dm = -f m, dp = -f p, dE = -f E) leaves u, v
// unchanged. Composition lanes are scaled (withdrawal) or incremented
// (injection) consistently with the element increments; the legacy total
// metal lane is updated when the composition sidecar is disabled.
// Magnetic face flux is never touched by a gas transaction (MHD owns it; the
// sink policy is recorded by the caller).
//
// Diagnostics use the actually committed increments, never requested ones.

#include <array>
#include <compare>
#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "cosmosim/physics/baryon_composition.hpp"

namespace cosmosim::core {
class SimulationState;
}

namespace cosmosim::physics {

enum class SourceChannel : std::uint8_t {
  kBhAccretionWithdrawal = 0,
  kBhSeedConversion = 1,
  kBhMerger = 2,
  kBhDynamicalFrictionReaction = 3,
  kAgnThermal = 4,
  kAgnKinetic = 5,
  kAgnWind = 6,
  kAgnJet = 7,
  kAgnMagnetic = 8,
  kAgnNozzle = 9,
  kRtMatterExchange = 10,
};

[[nodiscard]] const char* sourceChannelLabel(SourceChannel channel) noexcept;

struct SourceTransactionId {
  std::uint64_t source_id = 0U;        // stable BH particle ID / patch ID / event key
  std::uint64_t interval_index = 0U;   // accepted source-interval ordinal
  SourceChannel channel = SourceChannel::kBhAccretionWithdrawal;
  std::uint32_t sequence = 0U;         // deterministic per (source, interval, channel)
  auto operator<=>(const SourceTransactionId&) const = default;
};

// Integrated conserved increment in canonical code units (mass_code,
// momentum mass_code*velocity_code, energy mass_code*velocity_code^2).
struct ConservedIncrement {
  double mass_code = 0.0;
  std::array<double, 3> momentum_code{};
  std::array<double, 3> angular_momentum_code{};  // about the declared source point
  double total_energy_code = 0.0;                 // thermal + kinetic (no rest mass)
  std::array<double, k_chemical_element_count> element_mass_code{};
  double metal_mass_code = 0.0;                   // legacy total-metal lane

  ConservedIncrement& operator+=(const ConservedIncrement& rhs);
  [[nodiscard]] ConservedIncrement negated() const;
};

// Signed increment for one recipient (or donor) gas cell, addressed by
// stable identity. `local_row` is a prepared cache validated against
// gas_cell_id in prepare(); it is never trusted across epochs.
struct GasCellIncrement {
  std::uint64_t gas_cell_id = 0U;
  std::uint32_t local_row = 0U;
  ConservedIncrement delta{};
  // Withdrawal fraction (0..1) when the increment is a uniform removal; the
  // composition sidecar row is then scaled by (1 - fraction) and the
  // element_mass_code/metal increments are derived from the lanes in
  // prepare() (the caller sets mass/momentum/energy consistently).
  double withdrawal_fraction = 0.0;
};

enum class TransactionPhase : std::uint8_t {
  kPlanned = 0,
  kPrepared = 1,
  kValidated = 2,
  kAuthorized = 3,
  kCommitted = 4,
  kRejected = 5,
};

struct SourceTransactionPlan {
  SourceTransactionId id{};
  std::uint64_t topology_epoch = 0U;   // SimulationState::cellIndexGeneration() at planning
  std::uint32_t owner_rank = 0U;
  double evaluation_scale_factor = 1.0;
  double intended_time_code = 0.0;
  std::vector<GasCellIncrement> cells;
  // Increment received by the source object (BH / disk / pending budget /
  // radiation) so that sum(cells) + source_increment == external_increment.
  ConservedIncrement source_increment{};
  // Explicitly declared external budget (e.g. AGN energy drawn from the BH
  // engine ledger, radiation energy absorbed). Zero for closed transfers.
  ConservedIncrement external_increment{};
  // Tolerances (relative to the transferred magnitude, with absolute floors).
  double relative_tolerance = 1.0e-12;
  std::uint64_t reserved_scratch_bytes = 0U;
};

struct TransactionFloors {
  double mass_floor_code = 1.0e-20;
  double specific_internal_energy_floor_code = 0.0;
};

struct SourceTransactionReport {
  TransactionPhase phase = TransactionPhase::kPlanned;
  std::string failure;               // empty unless rejected
  std::uint64_t plans = 0U;
  std::uint64_t touched_rows = 0U;
  ConservedIncrement committed_to_gas{};      // actual summed gas increments
  ConservedIncrement committed_to_sources{};  // actual summed source increments
  double max_relative_closure_residual = 0.0;
};

// Deterministic batch of plans sharing one legal mutation epoch.
class SourceTransactionBatch {
 public:
  explicit SourceTransactionBatch(TransactionFloors floors = {});

  void add(SourceTransactionPlan plan);
  [[nodiscard]] std::size_t size() const noexcept { return m_plans.size(); }
  [[nodiscard]] std::span<const SourceTransactionPlan> plans() const noexcept { return m_plans; }
  [[nodiscard]] TransactionPhase phase() const noexcept { return m_phase; }

  // Pure checks; returns false (and sets the report failure) on rejection.
  [[nodiscard]] bool prepare(const core::SimulationState& state);
  [[nodiscard]] bool validate();
  // `participants_agree` is the caller-supplied unanimous vote (allreduce of
  // the local verdict). Passing false rejects the batch on every rank.
  [[nodiscard]] bool authorize(bool participants_agree);
  // Apply aggregated deltas exactly once; returns touched rows (sorted).
  std::vector<std::uint32_t> commit(core::SimulationState& state);

  [[nodiscard]] const SourceTransactionReport& report() const noexcept { return m_report; }

 private:
  struct RowAggregate {
    std::uint64_t gas_cell_id = 0U;
    std::uint32_t row = 0U;
    ConservedIncrement delta{};
    double retained_fraction = 1.0;  // product of (1 - withdrawal_fraction)
  };
  void reject(std::string reason);

  TransactionFloors m_floors;
  std::vector<SourceTransactionPlan> m_plans;
  std::vector<RowAggregate> m_rows;
  TransactionPhase m_phase = TransactionPhase::kPlanned;
  SourceTransactionReport m_report;
};

}  // namespace cosmosim::physics
