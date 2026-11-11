#include "cosmosim/amr/amr_geometry_coverage.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace cosmosim::amr {
namespace {

[[nodiscard]] std::string patchLabel(const AmrPatchGeometryRecord& record) {
  return "AMR patch " + std::to_string(record.patch_id) + " (rank " +
      std::to_string(record.owning_rank) + ", level " + std::to_string(record.level) + ")";
}

[[nodiscard]] std::string formatDouble(double value) {
  std::ostringstream stream;
  stream.precision(17);
  stream << value;
  return stream.str();
}

[[nodiscard]] bool cornerLess(const AmrSignedLatticeCorner& lhs, const AmrSignedLatticeCorner& rhs) noexcept {
  return std::tie(lhs.index[0], lhs.index[1], lhs.index[2]) <
      std::tie(rhs.index[0], rhs.index[1], rhs.index[2]);
}

[[nodiscard]] std::int64_t checkedAddWeight(std::int64_t lhs, std::int64_t rhs) {
  if ((rhs > 0 && lhs > std::numeric_limits<std::int64_t>::max() - rhs) ||
      (rhs < 0 && lhs < std::numeric_limits<std::int64_t>::min() - rhs)) {
    throw std::overflow_error("AMR coverage proof: signed corner weight overflow");
  }
  return lhs + rhs;
}

// Snaps a lattice coordinate to its integer plane. Returns false when the
// coordinate is non-finite or further than the snap tolerance from a plane.
[[nodiscard]] bool snapToLattice(double coordinate, double span, std::int64_t& index) noexcept {
  if (!std::isfinite(coordinate) ||
      std::abs(coordinate) > static_cast<double>(k_amr_coverage_max_lattice_cells_per_axis) + 1.0) {
    return false;
  }
  const double rounded = std::round(coordinate);
  if (std::abs(coordinate - rounded) > amrCoverageSnapTolerance(span)) {
    return false;
  }
  index = static_cast<std::int64_t>(rounded);
  return true;
}

}  // namespace

bool isLegacyUndeclaredAmrPatchGeometry(const AmrPatchGeometryRecord& record) noexcept {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (record.origin_comoving[axis] != 0.0 || record.extent_comoving[axis] != 0.0 ||
        record.cell_dims[axis] != 0U) {
      return false;
    }
  }
  return true;
}

std::string diagnoseAmrPatchGeometry(const AmrPatchGeometryRecord& record) {
  const std::string label = patchLabel(record);
  if (record.level < 0) {
    return label + " has an invalid negative refinement level";
  }
  if (record.cell_count == 0U) {
    return label + " is not populated and cannot be a tiling leaf";
  }
  std::uint64_t product = 1U;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const std::uint32_t dim = record.cell_dims[axis];
    if (dim == 0U) {
      return label + " axis " + std::to_string(axis) + ": cell dimension is 0";
    }
    if (product > std::numeric_limits<std::uint64_t>::max() / dim) {
      return label + " axis " + std::to_string(axis) + ": cell-dimension product overflows";
    }
    product *= dim;
  }
  if (product != record.cell_count) {
    return label + ": cell dimensions " + std::to_string(record.cell_dims[0]) + "x" +
        std::to_string(record.cell_dims[1]) + "x" + std::to_string(record.cell_dims[2]) +
        " do not match its populated cell count " + std::to_string(record.cell_count);
  }
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const double origin = record.origin_comoving[axis];
    const double extent = record.extent_comoving[axis];
    const std::string where = label + " axis " + std::to_string(axis) + ": ";
    if (!std::isfinite(origin)) {
      return where + "non-finite origin " + formatDouble(origin);
    }
    if (!std::isfinite(extent) || !(extent > 0.0)) {
      return where + "extent " + formatDouble(extent) + " is not finite and positive";
    }
    const double upper = origin + extent;
    if (!std::isfinite(upper) || !(upper > origin)) {
      return where + "upper bound " + formatDouble(upper) + " (origin " + formatDouble(origin) +
          " + extent " + formatDouble(extent) + ") is non-finite or not above the origin";
    }
    const double width = extent / static_cast<double>(record.cell_dims[axis]);
    if (!std::isfinite(width) || !(width > 0.0)) {
      return where + "cell width " + formatDouble(width) + " is not finite and positive";
    }
  }
  return {};
}

double amrCoverageSnapTolerance(double span_lattice_units) noexcept {
  const double span = std::isfinite(span_lattice_units) ? std::abs(span_lattice_units) : 0.0;
  return 1.0e-6 + 64.0 * std::numeric_limits<double>::epsilon() * span;
}

AmrCoverageLattice makeAmrCoverageLattice(
    const std::array<double, 3>& lower_comoving,
    const std::array<double, 3>& upper_comoving,
    const std::array<double, 3>& finest_cell_width_comoving) {
  AmrCoverageLattice lattice;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const double lower = lower_comoving[axis];
    const double upper = upper_comoving[axis];
    const double unit = finest_cell_width_comoving[axis];
    const std::string where = "AMR coverage lattice axis " + std::to_string(axis) + ": ";
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(upper > lower)) {
      throw std::runtime_error(where + "degenerate or non-finite bounding box [" + formatDouble(lower) +
                               ", " + formatDouble(upper) + ")");
    }
    if (!std::isfinite(unit) || !(unit > 0.0)) {
      throw std::runtime_error(where + "finest cell width " + formatDouble(unit) +
                               " is not finite and positive");
    }
    const double span = (upper - lower) / unit;
    if (!std::isfinite(span) || span > static_cast<double>(k_amr_coverage_max_lattice_cells_per_axis)) {
      throw std::runtime_error(where + "bounding box spans " + formatDouble(span) +
                               " finest cells, above the admitted 2^40 lattice span");
    }
    std::int64_t cells = 0;
    if (!snapToLattice(span, span, cells) || cells < 1) {
      throw std::runtime_error(where + "bounding box extent " + formatDouble(upper - lower) +
                               " is not an integer multiple of the finest cell width " + formatDouble(unit) +
                               " (" + formatDouble(span) + " units): mixed-level geometry is not dyadically aligned");
    }
    lattice.lower_comoving[axis] = lower;
    lattice.unit_comoving[axis] = unit;
    lattice.cells[axis] = cells;
  }
  return lattice;
}

void appendAmrPatchSignedCorners(
    const AmrCoverageLattice& lattice,
    const AmrPatchGeometryRecord& record,
    std::vector<AmrSignedLatticeCorner>& corners) {
  std::array<std::int64_t, 3> lo{0, 0, 0};
  std::array<std::int64_t, 3> hi{0, 0, 0};
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    const double unit = lattice.unit_comoving[axis];
    const double span = static_cast<double>(lattice.cells[axis]);
    const double origin = record.origin_comoving[axis];
    const double upper = origin + record.extent_comoving[axis];
    const double q_lo = (origin - lattice.lower_comoving[axis]) / unit;
    const double q_hi = (upper - lattice.lower_comoving[axis]) / unit;
    const std::string where = patchLabel(record) + " axis " + std::to_string(axis) + ": ";
    if (!snapToLattice(q_lo, span, lo[axis]) || !snapToLattice(q_hi, span, hi[axis])) {
      throw std::runtime_error(
          where + "face coordinates [" + formatDouble(origin) + ", " + formatDouble(upper) +
          ") are not on the finest-cell lattice (lattice units " + formatDouble(q_lo) + ", " +
          formatDouble(q_hi) + "); the hierarchy is not index-aligned");
    }
    if (lo[axis] < 0 || hi[axis] > lattice.cells[axis] || !(hi[axis] > lo[axis])) {
      throw std::runtime_error(where + "lattice interval [" + std::to_string(lo[axis]) + ", " +
                               std::to_string(hi[axis]) + ") leaves the bounding box [0, " +
                               std::to_string(lattice.cells[axis]) + ")");
    }
    const std::int64_t dim = static_cast<std::int64_t>(record.cell_dims[axis]);
    if ((hi[axis] - lo[axis]) % dim != 0) {
      throw std::runtime_error(where + "cell width is not an integer number of finest-cell lattice units (" +
                               std::to_string(hi[axis] - lo[axis]) + " units over " + std::to_string(dim) +
                               " cells)");
    }
  }
  corners.reserve(corners.size() + 8U);
  for (unsigned mask = 0U; mask < 8U; ++mask) {
    AmrSignedLatticeCorner corner;
    unsigned upper_axes = 0U;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const bool upper = (mask & (1U << axis)) != 0U;
      corner.index[axis] = upper ? hi[axis] : lo[axis];
      upper_axes += upper ? 1U : 0U;
    }
    corner.weight = (upper_axes % 2U) == 0U ? 1 : -1;
    corners.push_back(corner);
  }
}

void aggregateAmrSignedCorners(std::vector<AmrSignedLatticeCorner>& corners) {
  std::sort(corners.begin(), corners.end(), cornerLess);
  std::size_t write = 0U;
  for (std::size_t read = 0U; read < corners.size();) {
    AmrSignedLatticeCorner merged = corners[read];
    std::size_t next = read + 1U;
    while (next < corners.size() && corners[next].index == merged.index) {
      merged.weight = checkedAddWeight(merged.weight, corners[next].weight);
      ++next;
    }
    if (merged.weight != 0) {
      corners[write++] = merged;
    }
    read = next;
  }
  corners.resize(write);
}

long double amrPatchLatticeVolume(const AmrCoverageLattice& lattice, const AmrPatchGeometryRecord& record) {
  long double volume = 1.0L;
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    volume *= static_cast<long double>(record.extent_comoving[axis]) /
        static_cast<long double>(lattice.unit_comoving[axis]);
  }
  return volume;
}

AmrTilingVerdict verifyAmrExactTiling(
    const AmrCoverageLattice& lattice,
    std::span<const AmrSignedLatticeCorner> aggregated,
    long double summed_lattice_volume) {
  std::array<AmrSignedLatticeCorner, 8> expected{};
  for (unsigned mask = 0U; mask < 8U; ++mask) {
    unsigned upper_axes = 0U;
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      const bool upper = (mask & (1U << axis)) != 0U;
      expected[mask].index[axis] = upper ? lattice.cells[axis] : 0;
      upper_axes += upper ? 1U : 0U;
    }
    expected[mask].weight = (upper_axes % 2U) == 0U ? 1 : -1;
  }
  std::sort(expected.begin(), expected.end(), cornerLess);
  if (!std::is_sorted(aggregated.begin(), aggregated.end(), cornerLess)) {
    throw std::logic_error("verifyAmrExactTiling requires an aggregated (sorted) corner measure");
  }

  struct Residual {
    std::array<std::int64_t, 3> index;
    std::int64_t observed;
    std::int64_t expected;
  };
  std::vector<Residual> samples;
  std::uint64_t residual = 0U;
  const auto note = [&](const std::array<std::int64_t, 3>& index, std::int64_t observed, std::int64_t wanted) {
    ++residual;
    if (samples.size() < 4U) samples.push_back(Residual{index, observed, wanted});
  };
  std::size_t a = 0U;
  std::size_t e = 0U;
  while (a < aggregated.size() || e < expected.size()) {
    if (e == expected.size() || (a < aggregated.size() && cornerLess(aggregated[a], expected[e]))) {
      note(aggregated[a].index, aggregated[a].weight, 0);
      ++a;
    } else if (a == aggregated.size() || cornerLess(expected[e], aggregated[a])) {
      note(expected[e].index, 0, expected[e].weight);
      ++e;
    } else {
      if (aggregated[a].weight != expected[e].weight) {
        note(aggregated[a].index, aggregated[a].weight, expected[e].weight);
      }
      ++a;
      ++e;
    }
  }

  AmrTilingVerdict verdict;
  verdict.residual_points = residual;
  verdict.exact_tiling = residual == 0U;
  if (verdict.exact_tiling) {
    return verdict;
  }
  const long double box_volume = static_cast<long double>(lattice.cells[0]) *
      static_cast<long double>(lattice.cells[1]) * static_cast<long double>(lattice.cells[2]);
  const long double relative = (summed_lattice_volume - box_volume) / box_volume;
  std::string classification;
  if (relative > 1.0e-12L) {
    classification = "populated leaves overlap (summed volume exceeds the bounding box)";
  } else if (relative < -1.0e-12L) {
    classification = "populated leaves leave a hole (summed volume below the bounding box)";
  } else {
    classification = "populated leaves overlap AND leave a hole of compensating volume "
                     "(summed volume equals the bounding box; a volume check alone would accept this)";
  }
  std::string text = "AMR leaf tiling is not an exact cover of its bounding box: " + classification + "; " +
      std::to_string(residual) + " lattice point(s) carry a coverage-measure residual, e.g.";
  for (const Residual& sample : samples) {
    text += " (" + formatDouble(lattice.lower_comoving[0] + static_cast<double>(sample.index[0]) * lattice.unit_comoving[0]) +
        ", " + formatDouble(lattice.lower_comoving[1] + static_cast<double>(sample.index[1]) * lattice.unit_comoving[1]) +
        ", " + formatDouble(lattice.lower_comoving[2] + static_cast<double>(sample.index[2]) * lattice.unit_comoving[2]) +
        ") net " + std::to_string(sample.observed) + " expected " + std::to_string(sample.expected) + ";";
  }
  text += " production AMR hydro supports only a rectangular, hole-free, non-overlapping leaf tiling";
  verdict.diagnostic = std::move(text);
  return verdict;
}

AmrTilingVerdict verifyAmrPatchSetTiling(std::span<const AmrPatchGeometryRecord> records) {
  AmrTilingVerdict verdict;
  if (records.empty()) {
    verdict.diagnostic = "AMR coverage proof received no populated patch";
    return verdict;
  }
  const double inf = std::numeric_limits<double>::infinity();
  std::array<double, 3> lower{inf, inf, inf};
  std::array<double, 3> upper{-inf, -inf, -inf};
  std::array<double, 3> finest{inf, inf, inf};
  for (const AmrPatchGeometryRecord& record : records) {
    const std::string issue = diagnoseAmrPatchGeometry(record);
    if (!issue.empty()) {
      verdict.diagnostic = issue;
      return verdict;
    }
    for (std::size_t axis = 0; axis < 3U; ++axis) {
      lower[axis] = std::min(lower[axis], record.origin_comoving[axis]);
      upper[axis] = std::max(upper[axis], record.origin_comoving[axis] + record.extent_comoving[axis]);
      finest[axis] = std::min(finest[axis], record.extent_comoving[axis] / static_cast<double>(record.cell_dims[axis]));
    }
  }
  try {
    const AmrCoverageLattice lattice = makeAmrCoverageLattice(lower, upper, finest);
    std::vector<AmrSignedLatticeCorner> corners;
    corners.reserve(records.size() * 8U);
    long double volume = 0.0L;
    for (const AmrPatchGeometryRecord& record : records) {
      appendAmrPatchSignedCorners(lattice, record, corners);
      volume += amrPatchLatticeVolume(lattice, record);
    }
    aggregateAmrSignedCorners(corners);
    return verifyAmrExactTiling(lattice, corners, volume);
  } catch (const std::exception& error) {
    verdict.diagnostic = error.what();
    return verdict;
  }
}

}  // namespace cosmosim::amr
