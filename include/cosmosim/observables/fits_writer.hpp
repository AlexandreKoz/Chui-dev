#pragma once

// CHUI-EXT-PHYS-REF-v1 O01/RA05: minimal FITS (NASA/IAU FITS standard 4.0)
// primary-HDU writer for native sky products: BITPIX = -64 (IEEE double,
// big-endian), NAXIS <= 4 with linear WCS per axis (CTYPEn, CUNITn, CRPIXn,
// CRVALn, CDELTn), BUNIT, and provenance as HISTORY cards (68 characters per
// card, long text split, never truncated silently). Header and data are each
// padded to 2880-byte blocks. The file is written to `<path>.partial` and
// renamed on success, so readers never see a half-written product. No
// external process is invoked (no shell).

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace cosmosim::observables {

struct FitsAxis {
  std::size_t length = 0;
  std::string ctype;   // e.g. "VRAD", "FREQ", "OFFSET_X" (max 8 characters per FITS keyword value is not required)
  std::string cunit;   // e.g. "m/s", "Hz", "Mpc", "rad"
  double crpix = 1.0;  // 1-based reference pixel
  double crval = 0.0;
  double cdelt = 1.0;
};

struct FitsImage {
  std::vector<FitsAxis> axes;   // axis 1 varies fastest in `data`
  std::vector<double> data;
  std::string bunit;
  std::vector<std::pair<std::string, std::string>> string_keys;  // extra 8-character keywords
  std::vector<std::pair<std::string, double>> real_keys;
  std::vector<std::string> history;
};

// Builds the header card image (multiple of 2880 bytes); exposed for tests.
[[nodiscard]] std::string fitsHeader(const FitsImage& image);
void writeFitsImage(const std::string& path, const FitsImage& image);

}  // namespace cosmosim::observables
