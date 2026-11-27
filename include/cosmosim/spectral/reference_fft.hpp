#pragma once

// CHUI-EXT-PHYS-REF-v1 reference transforms for post-processing observables
// and species statistics (N05, L01, C02). Serial, deterministic, power-of-two
// iterative radix-2 Cooley-Tukey; NOT a replacement for the distributed PM
// FFT backends. Forward convention X_k = sum_n x_n exp(-2 pi i k n / N);
// inverse includes 1/N. Twiddles are computed per stage with std::polar.

#include <complex>
#include <cstddef>
#include <vector>

namespace cosmosim::spectral {

[[nodiscard]] bool isPowerOfTwo(std::size_t n) noexcept;
// In-place 1-D transform of a contiguous strided view (stride in elements).
void fftInPlace(std::complex<double>* data, std::size_t n, std::size_t stride, bool inverse);
void fft1d(std::vector<std::complex<double>>& data, bool inverse);
// Row-major cube [ix][iy][iz] of side n (power of two).
void fft3d(std::vector<std::complex<double>>& data, std::size_t n, bool inverse);
// Row-major square [ix][iy] of side n (power of two).
void fft2d(std::vector<std::complex<double>>& data, std::size_t n, bool inverse);
// Signed integer frequency of index i on an n-point periodic grid.
[[nodiscard]] long fftFrequencyIndex(std::size_t i, std::size_t n) noexcept;

}  // namespace cosmosim::spectral
