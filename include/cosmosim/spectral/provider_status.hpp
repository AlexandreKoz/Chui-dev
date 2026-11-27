#pragma once

// CHUI-EXT-PHYS-REF-v1 S02: provider availability contract. A missing
// physical provider must never return a plausible-looking zero, so every
// read-only physics provider in the extended modules reports one of three
// states. Consumers must branch on `status`; `value` is meaningful only for
// kAvailable (and is exactly 0 for kPhysicalZero, e.g. vacuum or a species
// genuinely absent from the state).

#include <stdexcept>
#include <string>

namespace cosmosim::spectral {

enum class ProviderStatus : unsigned char {
  kAvailable = 0,     // value computed from qualified-or-declared inputs
  kPhysicalZero = 1,  // the physical quantity is exactly zero (vacuum, absent species)
  kUnavailable = 2,   // no provider/data: the quantity is UNKNOWN, not zero
};

struct ProviderValue {
  double value = 0.0;
  ProviderStatus status = ProviderStatus::kUnavailable;
  std::string reason;  // non-empty for kUnavailable

  [[nodiscard]] static ProviderValue available(double v) { return {v, ProviderStatus::kAvailable, {}}; }
  [[nodiscard]] static ProviderValue physicalZero() { return {0.0, ProviderStatus::kPhysicalZero, {}}; }
  [[nodiscard]] static ProviderValue unavailable(std::string why) {
    return {0.0, ProviderStatus::kUnavailable, std::move(why)};
  }
  [[nodiscard]] bool known() const noexcept { return status != ProviderStatus::kUnavailable; }
  // Returns the value or throws with the recorded reason (fail closed).
  [[nodiscard]] double require(const char* what) const {
    if (status == ProviderStatus::kUnavailable) {
      throw std::runtime_error(std::string(what) + ": provider unavailable (" + reason + ")");
    }
    return value;
  }
};

}  // namespace cosmosim::spectral
