#!/usr/bin/env bash
# Full pipeline, from raw DEM to maps and animation.
set -euo pipefail
export SCRATCH="${SCRATCH:-/tmp/hawea}"
export DEM_DIR="$SCRATCH/dem"

python3 src/fetch_dem.py
python3 src/build_terrain.py
python3 src/build_bathymetry.py
python3 src/slope_screening.py
python3 src/impulse_wave.py          # validation against the VAW-211 worked example
python3 src/sensitivity.py
python3 src/seiche_modes.py

for S in S1_moderate S2_large S3_extreme; do
  python3 src/run_scenario.py "$S" 1500
done

python3 src/make_maps.py
python3 src/make_animation.py S2_large S3_extreme
echo "done - see outputs/figures and outputs/*.mp4"
