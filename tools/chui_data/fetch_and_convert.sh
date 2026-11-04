#!/usr/bin/env bash
# Reproducible acquisition + conversion of the CHUI external physics data
# (D5, 2026-11-04). Nothing here builds or runs CHUI.
#
# Usage: tools/chui_data/fetch_and_convert.sh <output_dir>
# Requires: curl, python3 (>= 3.9), h5py + numpy (WSS09 conversion only),
#           git (Starburst99 sample outputs from pyathena).
# Every SOURCE file is checked against the SHA-256 recorded on 2026-11-04;
# a mismatch aborts (an upstream file changed: re-review before using it).
# Derived files are written with the checksums printed by the converters;
# compare them with tools/chui_data/data_sources.md §"2026-11-04 acquisition".
set -euo pipefail
OUT=${1:?usage: fetch_and_convert.sh <output_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
RAW="$OUT/raw"
mkdir -p "$RAW" "$OUT"

fetch() {  # url file sha256
  local url=$1 file=$2 sum=$3
  if [ ! -f "$RAW/$file" ]; then curl -fsSL --retry 3 -o "$RAW/$file" "$url"; fi
  echo "$sum  $RAW/$file" | sha256sum -c --quiet - || { echo "checksum mismatch: $file" >&2; exit 2; }
}

NUPY=https://raw.githubusercontent.com/NuGrid/NuPyCEE/master
fetch $NUPY/yield_tables/agb_and_massive_stars_nugrid_MESAonly_fryer12delay.txt \
  agb_and_massive_stars_nugrid_MESAonly_fryer12delay.txt f2edf5ea29ed12d602022ae3605892b7cd0d4642eaf46945d813b2bdc9f002a1
fetch $NUPY/yield_tables/sn1a_i99_W7.txt sn1a_i99_W7.txt 352bb0b3d885c411a217f8242c7be792f5055d7444c2c802ed6601b11a33b39f
fetch $NUPY/LICENSE NuPyCEE_LICENSE 1ef5f600c55dc71d774028185104339473ddecebdaeccb7e513ddc2d9d49a3b5
CK=https://raw.githubusercontent.com/yipihey/ChemistryKernels.jl/main/data
fetch $CK/fg20_treecool_eff_rescaled_heating_rates_068.dat fg20_treecool_eff_rescaled_heating_rates_068.dat \
  4bfc5d9bd6f0541addd01256ac20b5e6bc8a7aa5fc8849b9cbdf34cb0882b313
fetch https://virgodb.cosma.dur.ac.uk/swift-webstorage/CoolingTables/EAGLE/coolingtables.tar.gz \
  eagle_coolingtables.tar.gz 228ec92871e8170a181ca5fb80a591b33daec4475a30b9482419f127670a1543
if [ ! -d "$RAW/pyathena" ]; then
  GIT_LFS_SKIP_SMUDGE=1 git clone -q --depth 1 --filter=blob:none --no-checkout \
    https://github.com/jeonggyukim/pyathena "$RAW/pyathena"
  git -C "$RAW/pyathena" checkout -q HEAD -- data/sb99/Z014_M1E6_GenevaV00_logdt_10Gyr LICENSE.md
fi
SB=$RAW/pyathena/data/sb99/Z014_M1E6_GenevaV00_logdt_10Gyr
echo "24aaa7b2f8bfca859a40cc11de7296c7fc5b34663b5d26c673c0be41a7e0397f  $SB/Z014_M1E6_GenevaV00_logdt_10Gyr.power1" | sha256sum -c --quiet -
echo "33f661447dcf788e82e5cd4b700d073ecd615452580928877c1eb87b66bac854  $SB/Z014_M1E6_GenevaV00_logdt_10Gyr.input1" | sha256sum -c --quiet -

# D1 yields
python3 -I "$HERE/convert_nugrid_yields.py" \
  --massive "$RAW/agb_and_massive_stars_nugrid_MESAonly_fryer12delay.txt" --ia "$RAW/sn1a_i99_W7.txt" \
  --output "$OUT/nugrid_fryer12delay_w7_chui_ssp_yield_table_v1.txt" \
  --citation "Pignatari et al. 2016 ApJS 225 24; Ritter et al. 2018 MNRAS 480 538; Iwamoto et al. 1999 ApJS 125 439; NuPyCEE (Cote et al. 2017)" \
  --license "NuPyCEE BSD-3-Clause (c) 2014-2016 NuGrid Team" --net-yield-sum-tolerance 0.02 --ia-untracked exclude
python3 -I "$HERE/validate_ssp_yield_table.py" "$OUT/nugrid_fryer12delay_w7_chui_ssp_yield_table_v1.txt" \
  --massive "$RAW/agb_and_massive_stars_nugrid_MESAonly_fryer12delay.txt" --ia "$RAW/sn1a_i99_W7.txt"
# D2 winds
python3 -I "$HERE/convert_sb99_winds.py" --power1 "$SB/Z014_M1E6_GenevaV00_logdt_10Gyr.power1" \
  --input1 "$SB/Z014_M1E6_GenevaV00_logdt_10Gyr.input1" --metallicity 0.014 \
  --citation "Leitherer et al. 1999 ApJS 123 3; Leitherer et al. 2014 ApJS 212 14 (Geneva v00 tracks, Ekstrom et al. 2012); run by J.-G. Kim (pyathena)" \
  --license "pyathena MIT License (c) 2019 J.-G. Kim & C.-G. Kim; Starburst99 output" \
  --output "$OUT/sb99_z014_genevav00_chui_ssp_wind_table_v1.txt"
# D3 UVB (validated in place) and WSS09 metal cooling
python3 -I "$HERE/validate_treecool.py" "$RAW/fg20_treecool_eff_rescaled_heating_rates_068.dat" \
  --expect-sha256 4bfc5d9bd6f0541addd01256ac20b5e6bc8a7aa5fc8849b9cbdf34cb0882b313
mkdir -p "$RAW/eagle" && tar xzf "$RAW/eagle_coolingtables.tar.gz" -C "$RAW/eagle"
python3 -I "$HERE/convert_wss09_metal_cooling.py" "$RAW"/eagle/coolingtables/z_[0-9].[0-9][0-9][0-9].hdf5 \
  --output "$OUT/wss09_hm01_z0-8.989_chui_metal_cooling_table_v1.txt" \
  --license "No licence stated by the authors (WSS08 page; SWIFT webstorage mirror); not redistributed" \
  --uv-background "HM01+CMB (Haardt & Madau 2001), CLOUDY 07.02, optically thin, ionization equilibrium"
python3 -I "$HERE/validate_metal_cooling_table.py" "$OUT/wss09_hm01_z0-8.989_chui_metal_cooling_table_v1.txt" \
  --sources "$RAW"/eagle/coolingtables/z_[0-9].[0-9][0-9][0-9].hdf5
cp "$RAW/fg20_treecool_eff_rescaled_heating_rates_068.dat" "$OUT/"
sha256sum "$OUT"/*.txt "$OUT"/*.dat
