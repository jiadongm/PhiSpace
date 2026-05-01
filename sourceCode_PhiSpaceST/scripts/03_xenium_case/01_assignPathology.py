#!/usr/bin/env python3
"""
assignPathology.py
Assign pathology annotations from GeoJSON to Xenium cells via point-in-polygon.

The GeoJSON annotations are in Xenium micron space (verified: GeoJSON x range
[59, 10783] matches cell x_centroid range [113, 10873]). The tissue strip is
~45mm in y, but the Xenium FOV only covers y=[92, 3608] (~3.5mm). Only the
Tumor polygon (Feature 0, y=[20, 35662]) overlaps the FOV. All other
annotations (Lymphoid, Blood Vessels, etc.) are at y > 14000, well outside
the FOV.

Smaller (more specific) polygons are checked first, so if a cell falls inside
both a Lymphoid aggregate and the broader Tumor polygon, it gets the Lymphoid
label.

Usage:
    python3 assignPathology.py

Output:
    data/Xenium_NSCLC/cell_pathology_annotations.csv
"""

import json
import gzip
import csv
from shapely.geometry import Point, shape

DATA_DIR = "data/Xenium_NSCLC"
GEOJSON_PATH = f"{DATA_DIR}/Xenium_V1_humanLung_Cancer_FFPE_annotation.geojson"
CELLS_PATH = f"{DATA_DIR}/Xenium_V1_humanLung_Cancer_FFPE_cells.csv.gz"
OUTPUT_PATH = f"{DATA_DIR}/cell_pathology_annotations.csv"


def main():
    # ---- 1. Load GeoJSON annotations ----
    with open(GEOJSON_PATH) as f:
        geo = json.load(f)

    annotations = []
    for feat in geo["features"]:
        props = feat.get("properties", {})
        classif = props.get("classification", {})
        name = classif.get("name", "Unknown") if isinstance(classif, dict) else "Unknown"
        poly = shape(feat["geometry"])
        annotations.append((name, poly, poly.area))
        print(f"  Loaded: {name:30s}  area={poly.area:>12.0f}  bounds={poly.bounds}")

    # Sort by area (smallest first) so specific annotations override broad ones
    annotations.sort(key=lambda x: x[2])
    print(f"\nPriority order (smallest/most-specific first):")
    for name, _, area in annotations:
        print(f"  {name}: area={area:.0f}")

    # ---- 2. Load cell coordinates ----
    cells = []
    with gzip.open(CELLS_PATH, "rt") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cells.append((row["cell_id"], float(row["x_centroid"]), float(row["y_centroid"])))

    print(f"\nLoaded {len(cells)} cells")
    print(f"  x_centroid: [{min(c[1] for c in cells):.1f}, {max(c[1] for c in cells):.1f}]")
    print(f"  y_centroid: [{min(c[2] for c in cells):.1f}, {max(c[2] for c in cells):.1f}]")

    # ---- 3. Point-in-polygon assignment ----
    # GeoJSON coordinates are in Xenium micron space (no transformation needed)
    results = []
    counts = {}
    for i, (cid, x, y) in enumerate(cells):
        pt = Point(x, y)
        label = "Unannotated"
        for name, poly, _ in annotations:
            if poly.contains(pt):
                label = name
                break
        results.append((cid, label))
        counts[label] = counts.get(label, 0) + 1

        if (i + 1) % 50000 == 0:
            print(f"  Processed {i + 1}/{len(cells)} cells...")

    # ---- 4. Report ----
    print(f"\nPathology assignment:")
    for label, cnt in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {label}: {cnt} ({100 * cnt / len(cells):.1f}%)")

    # ---- 5. Save ----
    with open(OUTPUT_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cell_id", "pathology"])
        for cid, label in results:
            writer.writerow([cid, label])

    print(f"\nSaved to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
