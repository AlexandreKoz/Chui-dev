#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "cosmosim/core/simulation_state.hpp"

namespace cosmosim::core {

// V02 canonical restart/rollback state. It coordinates measured chemistry
// demand with optional AMR resource selection; physical fields and mandatory
// refinement remain owned by their existing sidecars/hierarchy.
inline constexpr std::string_view k_v02_control_module =
    "chui_v02_joint_error_resource_control_v1";
inline constexpr std::uint32_t k_v02_control_schema = 1U;

struct V02ControlState {
  std::uint64_t epoch = 0U;
  double chemistry_promotion_bias = 0.0;
  double last_shadow_species_error = 0.0;
  std::uint64_t last_shadow_samples = 0U;
  std::uint64_t chemistry_detailed_rows = 0U;
  std::uint64_t chemistry_total_rows = 0U;
  std::uint64_t chemistry_promoted_last_interval = 0U;
  std::uint64_t last_headroom_bytes = 0U;
  std::uint64_t chemistry_growth_reserve_bytes = 0U;
  std::uint64_t optional_amr_budget_bytes = 0U;
  std::uint64_t optional_amr_requests = 0U;
  std::uint64_t optional_amr_admitted = 0U;
  // Normalized optional-indicator severity. 1 is the configured threshold.
  double optional_benefit_floor = 1.0;
};

namespace adaptive_resource_control_detail {

[[nodiscard]] inline std::string field(
    std::string_view text, std::string_view name, bool required) {
  const std::string key = std::string(name) + "=";
  const std::size_t begin = text.find(key);
  if (begin == std::string_view::npos) {
    if (!required) return {};
    throw std::runtime_error(
        "V02 controller sidecar is missing field '" + std::string(name) + "'");
  }
  const std::size_t value_begin = begin + key.size();
  const std::size_t end = text.find('\n', value_begin);
  return std::string(text.substr(value_begin, end - value_begin));
}

[[nodiscard]] inline std::uint64_t u64(
    std::string_view text, std::string_view name, std::uint64_t fallback) {
  const std::string value = field(text, name, false);
  if (value.empty()) return fallback;
  if (value.front() == '-') {
    throw std::runtime_error(
        "V02 controller sidecar has a negative integer field '" +
        std::string(name) + "'");
  }
  std::size_t consumed = 0U;
  const auto parsed = std::stoull(value, &consumed);
  if (consumed != value.size()) {
    throw std::runtime_error(
        "V02 controller sidecar has an invalid integer field '" + std::string(name) + "'");
  }
  return parsed;
}

[[nodiscard]] inline double real(
    std::string_view text, std::string_view name, double fallback) {
  const std::string value = field(text, name, false);
  if (value.empty()) return fallback;
  std::size_t consumed = 0U;
  const double parsed = std::stod(value, &consumed);
  if (consumed != value.size() || !std::isfinite(parsed)) {
    throw std::runtime_error(
        "V02 controller sidecar has an invalid floating field '" + std::string(name) + "'");
  }
  return parsed;
}

}  // namespace adaptive_resource_control_detail

[[nodiscard]] inline V02ControlState loadV02ControlState(
    const SimulationState& state) {
  V02ControlState result;
  const ModuleSidecarBlock* block = state.sidecars.find(k_v02_control_module);
  if (block == nullptr) return result;
  if (block->schema_version != k_v02_control_schema) {
    throw std::runtime_error("unsupported V02 controller sidecar schema");
  }
  if (block->payload.empty()) {
    throw std::runtime_error("V02 controller sidecar has an empty payload");
  }
  const std::string text(
      reinterpret_cast<const char*>(block->payload.data()), block->payload.size());
  using namespace adaptive_resource_control_detail;
  if (field(text, "module", true) !=
      "chui_v02_joint_error_resource_control") {
    throw std::runtime_error("V02 controller sidecar module marker mismatch");
  }
  result.epoch = u64(text, "epoch", 0U);
  result.chemistry_promotion_bias = real(text, "chemistry_promotion_bias", 0.0);
  result.last_shadow_species_error = real(text, "last_shadow_species_error", 0.0);
  result.last_shadow_samples = u64(text, "last_shadow_samples", 0U);
  result.chemistry_detailed_rows = u64(text, "chemistry_detailed_rows", 0U);
  result.chemistry_total_rows = u64(text, "chemistry_total_rows", 0U);
  result.chemistry_promoted_last_interval =
      u64(text, "chemistry_promoted_last_interval", 0U);
  result.last_headroom_bytes = u64(text, "last_headroom_bytes", 0U);
  result.chemistry_growth_reserve_bytes =
      u64(text, "chemistry_growth_reserve_bytes", 0U);
  result.optional_amr_budget_bytes = u64(text, "optional_amr_budget_bytes", 0U);
  result.optional_amr_requests = u64(text, "optional_amr_requests", 0U);
  result.optional_amr_admitted = u64(text, "optional_amr_admitted", 0U);
  result.optional_benefit_floor = real(text, "optional_benefit_floor", 1.0);
  if (result.chemistry_promotion_bias < 0.0 ||
      result.last_shadow_species_error < 0.0 ||
      result.optional_benefit_floor < 1.0 ||
      result.chemistry_detailed_rows > result.chemistry_total_rows ||
      result.optional_amr_admitted > result.optional_amr_requests) {
    throw std::runtime_error("V02 controller sidecar violates its state invariants");
  }
  return result;
}

inline void storeV02ControlState(
    SimulationState& state, const V02ControlState& control) {
  if (!std::isfinite(control.chemistry_promotion_bias) ||
      control.chemistry_promotion_bias < 0.0 ||
      !std::isfinite(control.last_shadow_species_error) ||
      control.last_shadow_species_error < 0.0 ||
      !std::isfinite(control.optional_benefit_floor) ||
      control.optional_benefit_floor < 1.0 ||
      control.chemistry_detailed_rows > control.chemistry_total_rows ||
      control.optional_amr_admitted > control.optional_amr_requests) {
    throw std::invalid_argument("invalid V02 controller state");
  }
  std::ostringstream text;
  text.precision(17);
  text << "module=chui_v02_joint_error_resource_control\n"
       << "epoch=" << control.epoch << '\n'
       << "chemistry_promotion_bias=" << control.chemistry_promotion_bias << '\n'
       << "last_shadow_species_error=" << control.last_shadow_species_error << '\n'
       << "last_shadow_samples=" << control.last_shadow_samples << '\n'
       << "chemistry_detailed_rows=" << control.chemistry_detailed_rows << '\n'
       << "chemistry_total_rows=" << control.chemistry_total_rows << '\n'
       << "chemistry_promoted_last_interval="
       << control.chemistry_promoted_last_interval << '\n'
       << "last_headroom_bytes=" << control.last_headroom_bytes << '\n'
       << "chemistry_growth_reserve_bytes="
       << control.chemistry_growth_reserve_bytes << '\n'
       << "optional_amr_budget_bytes=" << control.optional_amr_budget_bytes << '\n'
       << "optional_amr_requests=" << control.optional_amr_requests << '\n'
       << "optional_amr_admitted=" << control.optional_amr_admitted << '\n'
       << "optional_benefit_floor=" << control.optional_benefit_floor << '\n';
  const std::string payload = text.str();
  ModuleSidecarBlock block;
  block.module_name = std::string(k_v02_control_module);
  block.schema_version = k_v02_control_schema;
  block.payload.resize(payload.size());
  for (std::size_t i = 0; i < payload.size(); ++i) {
    block.payload[i] = static_cast<std::byte>(payload[i]);
  }
  state.sidecars.upsert(std::move(block));
}

}  // namespace cosmosim::core
