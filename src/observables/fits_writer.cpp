// CHUI-EXT-PHYS-REF-v1 O01 FITS writer (source only, NOT compiled in pass 01).

#include "cosmosim/observables/fits_writer.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace cosmosim::observables {
namespace {

constexpr std::size_t k_block = 2880U;

[[nodiscard]] std::string pad80(std::string card) {
  if (card.size() > 80U) throw std::invalid_argument("FITS card longer than 80 characters: " + card);
  card.resize(80U, ' ');
  return card;
}

[[nodiscard]] std::string keyword(const std::string& key) {
  if (key.empty() || key.size() > 8U) throw std::invalid_argument("FITS keyword must have 1-8 characters: " + key);
  for (const char c : key) {
    const bool ok = (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '-' || c == '_';
    if (!ok) throw std::invalid_argument("FITS keyword has an illegal character: " + key);
  }
  std::string k = key;
  k.resize(8U, ' ');
  return k;
}

[[nodiscard]] std::string logicalCard(const std::string& key, bool value) {
  std::string v(20U, ' ');
  v[19] = value ? 'T' : 'F';
  return pad80(keyword(key) + "= " + v);
}

[[nodiscard]] std::string integerCard(const std::string& key, long long value) {
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%20lld", value);
  return pad80(keyword(key) + "= " + buffer);
}

[[nodiscard]] std::string realCard(const std::string& key, double value) {
  if (!std::isfinite(value)) throw std::invalid_argument("FITS real keyword must be finite: " + key);
  char buffer[32];
  std::snprintf(buffer, sizeof(buffer), "%20.13E", value);
  return pad80(keyword(key) + "= " + buffer);
}

[[nodiscard]] std::string stringCard(const std::string& key, const std::string& value) {
  std::string escaped;
  for (const char c : value) {
    if (c < 32 || c > 126) throw std::invalid_argument("FITS string values must be printable ASCII");
    escaped += c;
    if (c == '\'') escaped += '\'';
  }
  if (escaped.size() > 68U) throw std::invalid_argument("FITS string value too long for one card: " + key);
  std::string quoted = "'" + escaped;
  if (quoted.size() < 9U) quoted.resize(9U, ' ');
  quoted += "'";
  return pad80(keyword(key) + "= " + quoted);
}

}  // namespace

std::string fitsHeader(const FitsImage& image) {
  if (image.axes.empty() || image.axes.size() > 4U) throw std::invalid_argument("FITS image needs 1-4 axes");
  std::size_t count = 1;
  for (const FitsAxis& a : image.axes) {
    if (a.length == 0U) throw std::invalid_argument("FITS axis length must be > 0");
    count *= a.length;
  }
  if (count != image.data.size()) throw std::invalid_argument("FITS data size does not match the axes");
  std::string header;
  header += logicalCard("SIMPLE", true);
  header += integerCard("BITPIX", -64);
  header += integerCard("NAXIS", static_cast<long long>(image.axes.size()));
  for (std::size_t i = 0; i < image.axes.size(); ++i) {
    header += integerCard("NAXIS" + std::to_string(i + 1U), static_cast<long long>(image.axes[i].length));
  }
  for (std::size_t i = 0; i < image.axes.size(); ++i) {
    const std::string n = std::to_string(i + 1U);
    const FitsAxis& a = image.axes[i];
    if (!a.ctype.empty()) header += stringCard("CTYPE" + n, a.ctype);
    if (!a.cunit.empty()) header += stringCard("CUNIT" + n, a.cunit);
    header += realCard("CRPIX" + n, a.crpix);
    header += realCard("CRVAL" + n, a.crval);
    header += realCard("CDELT" + n, a.cdelt);
  }
  if (!image.bunit.empty()) header += stringCard("BUNIT", image.bunit);
  for (const auto& [key, value] : image.string_keys) header += stringCard(key, value);
  for (const auto& [key, value] : image.real_keys) header += realCard(key, value);
  for (const std::string& text : image.history) {
    for (std::size_t pos = 0; pos < text.size(); pos += 72U) {
      std::string chunk = text.substr(pos, 72U);
      for (char& c : chunk) {
        if (c < 32 || c > 126) c = ' ';
      }
      header += pad80("HISTORY " + chunk);
    }
  }
  header += pad80("END");
  header.resize((header.size() + k_block - 1U) / k_block * k_block, ' ');
  return header;
}

void writeFitsImage(const std::string& path, const FitsImage& image) {
  const std::string header = fitsHeader(image);
  const std::string partial = path + ".partial";
  {
    std::ofstream out(partial, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("cannot open FITS output '" + partial + "'");
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    std::size_t written = 0;
    for (const double v : image.data) {
      std::uint64_t bits = std::bit_cast<std::uint64_t>(v);
      if constexpr (std::endian::native == std::endian::little) {
        bits = ((bits & 0x00000000000000FFULL) << 56U) | ((bits & 0x000000000000FF00ULL) << 40U) |
               ((bits & 0x0000000000FF0000ULL) << 24U) | ((bits & 0x00000000FF000000ULL) << 8U) |
               ((bits & 0x000000FF00000000ULL) >> 8U) | ((bits & 0x0000FF0000000000ULL) >> 24U) |
               ((bits & 0x00FF000000000000ULL) >> 40U) | ((bits & 0xFF00000000000000ULL) >> 56U);
      }
      char bytes[8];
      std::memcpy(bytes, &bits, 8U);
      out.write(bytes, 8);
      written += 8U;
    }
    const std::size_t pad = (k_block - written % k_block) % k_block;
    const std::string zeros(pad, '\0');
    out.write(zeros.data(), static_cast<std::streamsize>(zeros.size()));
    if (!out) throw std::runtime_error("FITS write failed for '" + partial + "'");
  }
  std::filesystem::rename(partial, path);
}

}  // namespace cosmosim::observables
