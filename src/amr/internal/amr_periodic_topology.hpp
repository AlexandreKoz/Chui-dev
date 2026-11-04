#pragma once

// CHUI Part 3 Checkpoint A — periodic AMR patch topology helpers.
//
// Box convention: [0, L_axis) per axis (same as core::canonicalPeriodicPosition
// and the TreePM domain). Every real AMR patch lies inside the box. Periodic
// neighbour relations are expressed WITHOUT materializing translated patch
// copies:
//
//  * face contact: the normal coordinates of two patch faces are compared
//    with the minimum-image difference (a lower face at 0 touches an upper
//    face at L);
//  * point / volume queries (ghost probes, fine-to-coarse ghost volumes) are
//    canonicalized into [0, L) before the ordinary containment test.
//
// An axis with extent <= 0 (or a frame that is disabled) is non-periodic, so
// isolated modes keep their exact historical geometry. These helpers are
// header-only, allocation-free and used by both the orchestrator and the
// ghost filler, so the two can never disagree about wrapped topology.

#include <algorithm>
#include <array>
#include <cmath>

namespace cosmosim::amr::internal {

struct AmrPeriodicFrame {
  bool enabled = false;
  std::array<double, 3> extent_comoving{0.0, 0.0, 0.0};

  [[nodiscard]] bool axisPeriodic(std::size_t axis) const noexcept {
    return enabled && extent_comoving[axis] > 0.0 && std::isfinite(extent_comoving[axis]);
  }
};

[[nodiscard]] inline AmrPeriodicFrame makeAmrPeriodicFrame(
    bool periodic, const std::array<double, 3>& extent_comoving) noexcept {
  AmrPeriodicFrame frame;
  frame.enabled = periodic;
  frame.extent_comoving = extent_comoving;
  return frame;
}

// Canonical coordinate in [0, L). A coordinate within a relative 1e-12 of L is
// mapped to 0 so that a probe placed exactly on the upper box face resolves to
// the first cell of the wrapped neighbour rather than escaping the box.
[[nodiscard]] inline double wrapCoordinate(double x, double length) noexcept {
  if (!(length > 0.0) || !std::isfinite(x)) {
    return x;
  }
  double wrapped = std::fmod(x, length);
  if (wrapped < 0.0) {
    wrapped += length;
  }
  if (wrapped >= length * (1.0 - 1.0e-12)) {
    wrapped = 0.0;
  }
  return wrapped;
}

[[nodiscard]] inline std::array<double, 3> wrapPoint(
    const AmrPeriodicFrame& frame, std::array<double, 3> point) noexcept {
  for (std::size_t axis = 0; axis < 3U; ++axis) {
    if (frame.axisPeriodic(axis)) {
      point[axis] = wrapCoordinate(point[axis], frame.extent_comoving[axis]);
    }
  }
  return point;
}

// Minimum-image difference a - b along one axis.
[[nodiscard]] inline double minImageDifference(
    const AmrPeriodicFrame& frame, std::size_t axis, double a, double b) noexcept {
  double delta = a - b;
  if (frame.axisPeriodic(axis)) {
    const double length = frame.extent_comoving[axis];
    delta -= length * std::round(delta / length);
  }
  return delta;
}

// Shift that moves the axis-aligned interval [lower, upper) into the box when
// it lies (entirely) outside it. Ghost volumes are exactly one target cell
// wide and start on a cell boundary that coincides with the box face, so a
// single shift by +-L is always sufficient.
[[nodiscard]] inline double intervalWrapShift(
    const AmrPeriodicFrame& frame, std::size_t axis, double lower, double upper) noexcept {
  if (!frame.axisPeriodic(axis)) {
    return 0.0;
  }
  const double length = frame.extent_comoving[axis];
  const double tolerance = 1.0e-12 * std::max(1.0, length);
  if (upper <= tolerance) {
    return length;
  }
  if (lower >= length - tolerance) {
    return -length;
  }
  return 0.0;
}

// True when the patch spans the whole periodic axis (in-patch wrap applies).
[[nodiscard]] inline bool spansPeriodicAxis(
    const AmrPeriodicFrame& frame, std::size_t axis, double extent, double tolerance) noexcept {
  return frame.axisPeriodic(axis) &&
      std::abs(extent - frame.extent_comoving[axis]) <= tolerance;
}

}  // namespace cosmosim::amr::internal
