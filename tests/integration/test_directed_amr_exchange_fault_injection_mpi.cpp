// CHUI patch 53 — P0-52-01 directed AMR exchange collective-safety source
// tests (NOT RUN — USER-RESERVED for the validation campaign).
//
// Run with 3 or 4 ranks. Slab decomposition along x, one 2x1x1 patch per rank
// except the last rank, which owns no patch (zero-data participant). Each
// scenario must complete on EVERY rank (no hang) and every rank must reach the
// same verdict (all succeed or all throw); the verdict is cross-checked with
// an MPI_Allreduce after each scenario. A hang is the failure mode under test,
// so the ctest TIMEOUT is the oracle for deadlock.
//
// Scenarios:
//  1. fault-free exchange (baseline; deterministic descriptor order);
//  2. consumer callback throws on rank 1 in the first streamed round;
//  3. producer callback throws on rank 0;
//  4. admission callback refuses on the highest patch-owning rank;
//  5. per-rank MemoryGovernor with a 64-byte hard limit on rank 0 only
//     (unilateral control-plane / descriptor admission refusal);
//  6. one-record transport rounds (multi-round streaming, int-safe counts);
//  7. periodic topology (rank 0 and the last patch rank meet across the wrap);
//  8. flux-register exchange with a stale-owner record injected on rank 1.
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {
namespace par = cosmosim::parallel;
namespace core = cosmosim::core;

#if COSMOSIM_ENABLE_MPI
int g_rank = 0;
int g_size = 1;

par::AmrPatchPayloadRecord makePatch(int rank, int patch_ranks) {
  par::AmrPatchPayloadRecord record;
  record.patch_id = 10000U + static_cast<std::uint64_t>(rank);
  record.morton_key = static_cast<std::uint64_t>(rank + 1);
  record.owner_rank = rank;
  record.cell_count = 2U;
  const double width = 1.0 / static_cast<double>(patch_ranks);
  record.origin_x_comoving = width * static_cast<double>(rank);
  record.extent_x_comoving = width;
  record.extent_y_comoving = 1.0;
  record.extent_z_comoving = 1.0;
  record.cell_dim_x = 2U;
  record.cell_dim_y = 1U;
  record.cell_dim_z = 1U;
  record.decomposition_epoch = 53U;
  record.cell_mass_sum_code = 2.0;
  record.gas_internal_energy_sum_code = 3.0;
  par::validateAmrPatchPayloadRecord(record);
  return record;
}

par::AmrPatchCellPayloadRecord makeCell(const par::AmrPatchPayloadRecord& patch, std::uint32_t offset) {
  par::AmrPatchCellPayloadRecord record;
  record.patch_id = patch.patch_id;
  record.owner_rank = patch.owner_rank;
  record.local_cell_offset = offset;
  record.center_x_comoving = patch.origin_x_comoving + patch.extent_x_comoving * (0.25 + 0.5 * offset);
  record.center_y_comoving = 0.5;
  record.center_z_comoving = 0.5;
  record.mass_code = 1.0;
  record.gas_cell_id = 70000U + 10U * static_cast<std::uint64_t>(patch.owner_rank) + offset;
  record.parent_particle_id = record.gas_cell_id + 1U;
  record.density_code = 1.0;
  record.pressure_code = 1.0;
  record.internal_energy_code = 1.0;
  record.temperature_code = 1.0;
  record.sound_speed_code = 1.0;
  par::validateAmrPatchCellPayloadRecord(record);
  return record;
}

struct Faults {
  int consumer_rank = -1;
  int producer_rank = -1;
  int admission_rank = -1;
  int tight_governor_rank = -1;
  std::size_t round_limit_bytes = 0U;
  bool periodic = false;
};

// Returns 1 when this rank's exchange threw, 0 otherwise.
int runDirected(const Faults& faults) {
  const int patch_ranks = g_size - 1;
  std::vector<par::AmrPatchPayloadRecord> local;
  if (g_rank < patch_ranks) local.push_back(makePatch(g_rank, patch_ranks));
  core::MemoryGovernorPolicy policy;
  if (g_rank == faults.tight_governor_rank) policy.hard_limit_bytes = 64U;
  core::MemoryGovernor governor(policy);
  std::vector<par::AmrPatchCellPayloadRecord> received;
  try {
    const par::MpiContext mpi;
    const auto exchange = par::executeBlockingDirectedAmrPatchPayloadExchange(
        mpi, local,
        [&](std::span<const par::AmrPatchBoundaryCellRequest> requests, std::uint64_t first, std::size_t max,
            std::vector<par::AmrPatchCellPayloadRecord>& out) {
          if (g_rank == faults.producer_rank) throw std::runtime_error("injected producer failure");
          std::vector<par::AmrPatchCellPayloadRecord> logical;
          // 2x1x1 patch, one ghost layer: the x-lower face cell is offset 0,
          // the x-upper face cell offset 1 (the only faces shared by slabs).
          for (const auto& request : requests) {
            const auto lower = static_cast<std::uint8_t>(par::AmrPatchBoundaryFace::kXLower);
            const auto upper = static_cast<std::uint8_t>(par::AmrPatchBoundaryFace::kXUpper);
            if ((request.boundary_face_mask & lower) != 0U) logical.push_back(makeCell(local.front(), 0U));
            if ((request.boundary_face_mask & upper) != 0U) logical.push_back(makeCell(local.front(), 1U));
          }
          for (std::size_t i = 0; i < max && first + i < logical.size(); ++i) out.push_back(logical[first + i]);
        },
        [&](int, std::span<const par::AmrPatchCellPayloadRecord> records) {
          if (g_rank == faults.consumer_rank) throw std::runtime_error("injected consumer failure");
          received.insert(received.end(), records.begin(), records.end());
        },
        [&](int, std::span<const par::AmrPatchPayloadRecord>, std::span<const par::AmrPatchBoundaryCellRequest>,
            std::uint64_t) {
          if (g_rank == faults.admission_rank) throw std::runtime_error("injected admission refusal");
        },
        faults.round_limit_bytes, 53U,
        par::AmrInterfaceTopologyOptions{
            .periodic_box_extent_comoving = faults.periodic ? std::array<double, 3>{1.0, 1.0, 1.0}
                                                            : std::array<double, 3>{0.0, 0.0, 0.0},
            .ghost_layers = 1U},
        &governor);
    // Deterministic order of the agreed descriptors (ascending peer rank).
    for (std::size_t i = 1; i < exchange.patch_payloads_received.size(); ++i) {
      assert(exchange.patch_payloads_received[i - 1].owner_rank <= exchange.patch_payloads_received[i].owner_rank);
    }
    return 0;
  } catch (const std::exception&) {
    return 1;
  }
}

void expectAgreement(int threw, int expect_throw, const char* label) {
  int sum = 0;
  MPI_Allreduce(&threw, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (!(sum == 0 || sum == g_size)) {
    std::fprintf(stderr, "[%s] ranks disagree on the exchange verdict (%d/%d threw)\n", label, sum, g_size);
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
  if ((sum == g_size) != (expect_throw != 0)) {
    std::fprintf(stderr, "[%s] unexpected collective verdict\n", label);
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
}

int runFlux(bool inject_stale) {
  std::vector<par::AmrFluxRegisterPayloadRecord> local;
  par::AmrFluxRegisterPayloadRecord record;
  record.register_key = 99000U + static_cast<std::uint64_t>(g_rank);
  record.coarse_patch_id = 10000U;
  record.coarse_gas_cell_id = 70000U;
  record.level = 1U;
  record.source_rank = g_rank;
  record.owner_rank = (g_rank + 1) % g_size;
  record.gas_cell_identity_generation = 1U;
  record.patch_geometry_generation = 1U;
  record.face_area_comov = 0.25;
  record.coarse_area_comov = 0.25;
  record.fine_area_comov = 0.25;
  record.dt_code = 0.001;
  record.coarse_face_count = 1U;
  record.fine_face_count = 1U;
  local.push_back(record);
  if (inject_stale && g_rank == 1) local.front().source_rank = 0;  // stale metadata on one rank
  try {
    const par::MpiContext mpi;
    const auto result = par::executeGovernedAmrFluxRegisterPayloadExchange(mpi, local, 53U, nullptr);
    (void)result;
    return 0;
  } catch (const std::exception&) {
    return 1;
  }
}
#endif

}  // namespace

int main() {
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &g_size);
  if (g_size < 3) {
    MPI_Finalize();
    return 0;
  }
  expectAgreement(runDirected(Faults{}), 0, "baseline");
  expectAgreement(runDirected(Faults{.consumer_rank = 1}), 1, "consumer failure");
  expectAgreement(runDirected(Faults{.producer_rank = 0}), 1, "producer failure");
  expectAgreement(runDirected(Faults{.admission_rank = g_size - 2}), 1, "admission refusal");
  expectAgreement(runDirected(Faults{.tight_governor_rank = 0}), 1, "governor refusal");
  expectAgreement(runDirected(Faults{.round_limit_bytes = sizeof(par::AmrPatchCellPayloadRecord)}), 0,
                  "one-record rounds");
  expectAgreement(runDirected(Faults{.periodic = true}), 0, "periodic topology");
  expectAgreement(runFlux(false), 0, "flux baseline");
  expectAgreement(runFlux(true), 1, "flux stale metadata");
  MPI_Finalize();
#endif
  return 0;
}
