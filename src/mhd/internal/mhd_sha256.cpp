#include "mhd_sha256.hpp"

#include <array>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace cosmosim::mhd::internal {
namespace {

constexpr std::array<std::uint32_t, 64> k_round_constants{
    0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U,
    0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
    0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U, 0xe49b69c1U, 0xefbe4786U,
    0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
    0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
    0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
    0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
    0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
    0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU,
    0x5b9cca4fU, 0x682e6ff3U, 0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
    0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

constexpr std::array<std::uint32_t, 8> k_initial_hash{
    0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
    0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U};

std::uint32_t rotateRight(std::uint32_t value, unsigned count) noexcept {
  return (value >> count) | (value << (32U - count));
}

void compressBlock(const std::uint8_t* block, std::array<std::uint32_t, 8>& hash) noexcept {
  std::array<std::uint32_t, 64> schedule{};
  for (std::size_t index = 0; index < 16U; ++index) {
    schedule[index] = (static_cast<std::uint32_t>(block[index * 4U]) << 24U) |
        (static_cast<std::uint32_t>(block[index * 4U + 1U]) << 16U) |
        (static_cast<std::uint32_t>(block[index * 4U + 2U]) << 8U) |
        static_cast<std::uint32_t>(block[index * 4U + 3U]);
  }
  for (std::size_t index = 16U; index < 64U; ++index) {
    const std::uint32_t s0 = rotateRight(schedule[index - 15U], 7U) ^
        rotateRight(schedule[index - 15U], 18U) ^ (schedule[index - 15U] >> 3U);
    const std::uint32_t s1 = rotateRight(schedule[index - 2U], 17U) ^
        rotateRight(schedule[index - 2U], 19U) ^ (schedule[index - 2U] >> 10U);
    schedule[index] = schedule[index - 16U] + s0 + schedule[index - 7U] + s1;
  }
  std::uint32_t a = hash[0];
  std::uint32_t b = hash[1];
  std::uint32_t c = hash[2];
  std::uint32_t d = hash[3];
  std::uint32_t e = hash[4];
  std::uint32_t f = hash[5];
  std::uint32_t g = hash[6];
  std::uint32_t h = hash[7];
  for (std::size_t index = 0; index < 64U; ++index) {
    const std::uint32_t sigma1 =
        rotateRight(e, 6U) ^ rotateRight(e, 11U) ^ rotateRight(e, 25U);
    const std::uint32_t choose = (e & f) ^ (~e & g);
    const std::uint32_t temp1 = h + sigma1 + choose + k_round_constants[index] + schedule[index];
    const std::uint32_t sigma0 =
        rotateRight(a, 2U) ^ rotateRight(a, 13U) ^ rotateRight(a, 22U);
    const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
    const std::uint32_t temp2 = sigma0 + majority;
    h = g;
    g = f;
    f = e;
    e = d + temp1;
    d = c;
    c = b;
    b = a;
    a = temp1 + temp2;
  }
  hash[0] += a;
  hash[1] += b;
  hash[2] += c;
  hash[3] += d;
  hash[4] += e;
  hash[5] += f;
  hash[6] += g;
  hash[7] += h;
}

std::string toHex(const std::array<std::uint32_t, 8>& hash) {
  static constexpr char k_digits[] = "0123456789abcdef";
  std::string out;
  out.reserve(64U);
  for (const std::uint32_t word : hash) {
    for (int shift = 28; shift >= 0; shift -= 4) {
      out.push_back(k_digits[(word >> static_cast<unsigned>(shift)) & 0xFU]);
    }
  }
  return out;
}

}  // namespace

std::string sha256Hex(std::string_view bytes) {
  std::array<std::uint32_t, 8> hash = k_initial_hash;
  const std::size_t length = bytes.size();
  const std::size_t full_blocks = length / 64U;
  for (std::size_t block = 0; block < full_blocks; ++block) {
    compressBlock(reinterpret_cast<const std::uint8_t*>(bytes.data()) + block * 64U, hash);
  }
  // Remainder plus padding.
  const std::size_t remainder = length % 64U;
  std::array<std::uint8_t, 128> tail{};
  for (std::size_t index = 0; index < remainder; ++index) {
    tail[index] = static_cast<std::uint8_t>(bytes[full_blocks * 64U + index]);
  }
  tail[remainder] = 0x80U;
  const std::size_t padded_length = remainder < 56U ? 64U : 128U;
  const std::uint64_t bit_length = static_cast<std::uint64_t>(length) * 8U;
  for (std::size_t index = 0; index < 8U; ++index) {
    tail[padded_length - 1U - index] =
        static_cast<std::uint8_t>((bit_length >> (8U * index)) & 0xFFU);
  }
  compressBlock(tail.data(), hash);
  if (padded_length == 128U) {
    compressBlock(tail.data() + 64U, hash);
  }
  return toHex(hash);
}

std::string sha256FileHex(const std::filesystem::path& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) {
    throw std::runtime_error("mhd sha256: cannot open " + path.string());
  }
  std::string contents((std::istreambuf_iterator<char>(stream)),
                       std::istreambuf_iterator<char>());
  if (!stream.eof() && stream.fail()) {
    throw std::runtime_error("mhd sha256: cannot read " + path.string());
  }
  return sha256Hex(contents);
}

}  // namespace cosmosim::mhd::internal
