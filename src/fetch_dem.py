"""
Download the Copernicus GLO-30 DEM tiles covering Lake Hawea.

Copernicus DEM (ESA / Airbus) is free to use and redistribute with attribution.
Tiles are public COGs on AWS Open Data - no credentials needed.
"""
import os, sys, urllib.request

DEST = os.environ.get("DEM_DIR", "/tmp/hawea_dem")
BASE = ("https://copernicus-dem-30m.s3.amazonaws.com/"
        "Copernicus_DSM_COG_10_{t}_DEM/Copernicus_DSM_COG_10_{t}_DEM.tif")
TILES = ["S45_00_E169_00", "S45_00_E168_00", "S44_00_E169_00"]

def main():
    os.makedirs(DEST, exist_ok=True)
    for t in TILES:
        out = os.path.join(DEST, f"cop30_{t}.tif")
        if os.path.exists(out) and os.path.getsize(out) > 1e6:
            print(f"have {out}"); continue
        url = BASE.format(t=t)
        print(f"fetching {t} ...", flush=True)
        urllib.request.urlretrieve(url, out)
        print(f"  -> {out} ({os.path.getsize(out)/1e6:.1f} MB)")
    print(f"\nSet SCRATCH so build_terrain.py finds them, e.g.:\n  export SCRATCH={os.path.dirname(DEST)}")

if __name__ == "__main__":
    main()
