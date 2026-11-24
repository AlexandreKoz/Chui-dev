#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

namespace cosmosim::mhd::internal {

// SHA-256 of a byte string, lowercase hexadecimal.
[[nodiscard]] std::string sha256Hex(std::string_view bytes);

// SHA-256 of a file's contents, lowercase hexadecimal. Throws std::runtime_error
// when the file cannot be opened or read completely.
[[nodiscard]] std::string sha256FileHex(const std::filesystem::path& path);

}  // namespace cosmosim::mhd::internal
