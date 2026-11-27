// CHUI-EXT-PHYS-REF-v1 reference radix-2 FFT (source only, NOT compiled in pass 01).

#include "cosmosim/spectral/reference_fft.hpp"

#include <stdexcept>
#include <utility>

namespace cosmosim::spectral {

bool isPowerOfTwo(std::size_t n) noexcept { return n != 0U && (n & (n - 1U)) == 0U; }

void fftInPlace(std::complex<double>* data, std::size_t n, std::size_t stride, bool inverse) {
  if (!isPowerOfTwo(n)) throw std::invalid_argument("reference FFT length must be a power of two");
  if (n == 1U) return;
  // Bit-reversal permutation.
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1U;
    for (; (j & bit) != 0U; bit >>= 1U) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(data[i * stride], data[j * stride]);
  }
  const double pi = 3.141592653589793238462643383279502884;
  for (std::size_t len = 2; len <= n; len <<= 1U) {
    const double angle = (inverse ? 2.0 : -2.0) * pi / static_cast<double>(len);
    for (std::size_t start = 0; start < n; start += len) {
      for (std::size_t m = 0; m < len / 2U; ++m) {
        const std::complex<double> w = std::polar(1.0, angle * static_cast<double>(m));
        std::complex<double>& lo = data[(start + m) * stride];
        std::complex<double>& hi = data[(start + m + len / 2U) * stride];
        const std::complex<double> t = w * hi;
        hi = lo - t;
        lo = lo + t;
      }
    }
  }
  if (inverse) {
    const double scale = 1.0 / static_cast<double>(n);
    for (std::size_t i = 0; i < n; ++i) data[i * stride] *= scale;
  }
}

void fft1d(std::vector<std::complex<double>>& data, bool inverse) { fftInPlace(data.data(), data.size(), 1U, inverse); }

void fft2d(std::vector<std::complex<double>>& data, std::size_t n, bool inverse) {
  if (data.size() != n * n) throw std::invalid_argument("fft2d size mismatch");
  for (std::size_t ix = 0; ix < n; ++ix) fftInPlace(data.data() + ix * n, n, 1U, inverse);
  for (std::size_t iy = 0; iy < n; ++iy) fftInPlace(data.data() + iy, n, n, inverse);
}

void fft3d(std::vector<std::complex<double>>& data, std::size_t n, bool inverse) {
  if (data.size() != n * n * n) throw std::invalid_argument("fft3d size mismatch");
  for (std::size_t ix = 0; ix < n; ++ix) {
    for (std::size_t iy = 0; iy < n; ++iy) fftInPlace(data.data() + (ix * n + iy) * n, n, 1U, inverse);
  }
  for (std::size_t ix = 0; ix < n; ++ix) {
    for (std::size_t iz = 0; iz < n; ++iz) fftInPlace(data.data() + ix * n * n + iz, n, n, inverse);
  }
  for (std::size_t iy = 0; iy < n; ++iy) {
    for (std::size_t iz = 0; iz < n; ++iz) fftInPlace(data.data() + iy * n + iz, n, n * n, inverse);
  }
}

long fftFrequencyIndex(std::size_t i, std::size_t n) noexcept {
  return i <= n / 2U ? static_cast<long>(i) : static_cast<long>(i) - static_cast<long>(n);
}

}  // namespace cosmosim::spectral
