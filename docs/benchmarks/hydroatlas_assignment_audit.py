#!/usr/bin/env python
"""Inputs for the HydroATLAS containing-subbasin rule (user decisions 2026-11-19) and the per-gage topology audit.

The rule (applied in R by `apply_containing_subbasin_rule()`, R/aggregate_hydroatlas_metadata.R):
    a gage gets its containing level-12 subbasin c ALONE when c is not a HydroBASINS headwater
    (some subbasin has NEXT_DOWN = c) AND SUB_AREA(c) > 1.5 x the official polygon area;
    every other gage keeps c plus everything upstream.
The official polygon area is `geom_area_km2` of the agency outline (watershed_geom_source gagesii /
wsc_eccc) in the deposited unsimplified boundary layer. Gages with no agency outline (the HydroBASINS
fallbacks and 05KH009) never trigger the rule.

This tool writes three files to --out:
  basinatlas_attributes.parquet          every BasinATLAS v10 level-12 attribute without geometry
                                         (HYBAS_ID / NEXT_DOWN as strings), the attribute source for the R step;
  hydroatlas_assignment_audit.csv        one row per attribute-table gage (8,014): the official area, c, the rule
                                         decision, and the topology check against the agency polygon, before
                                         and after the rule, with `assignment_flag`;
  hydroatlas_assignment_audit_provenance.json   input hashes and the headline counts.

Topology check (acceptance standard, user decisions 2026-11-19; never IoU). All areas in ESRI:102008.
For an assigned set S (c + upstream, or c alone) and the agency polygon P:
    exc_up  = area of S outside P, not counting c itself (resolution excess inside c is accepted)
    deficit = share of P not covered by S
    defect  = exc_up > 0.5 P, or deficit > 0.5 unless the gage is smaller than c (P < SUB_AREA(c)) and the
              missed area does not lie mostly upstream of c (a small watershed straddling a divide is accepted).
assignment_flag = no_official_polygon | topology_defect | ok.
Agency polygons with more than 200,000 vertices (16) are simplified at 100 m before the overlay, as in the
2026-11-19 review (docs/plans/2026-11-19-hydroatlas-outlet-review/).

Environment: pandas, pyarrow, numpy, shapely >= 2, pyproj, pyogrio (the repo .venv lacks shapely):
    uv venv --python 3.12 /tmp/ha_env && uv pip install --python /tmp/ha_env/bin/python pandas pyarrow numpy shapely pyproj pyogrio
Usage:
    python docs/benchmarks/hydroatlas_assignment_audit.py --out ~/HISSS_data/hydroatlas-assignment-fix-19nov2026/inputs   # ~40 min
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import sqlite3
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyogrio
import pyproj
import shapely

HOME = os.path.expanduser("~")
R4 = f"{HOME}/Downloads/Signatures_v3_28oct2026/resource4_geometry_attributes"
DEF_BASINATLAS = f"{HOME}/Downloads/geometry_rebuild/basinAt_NorAm_polys.gpkg"
DEF_ATTR = f"{R4}/hisss_watershed_hydroatlas_attributes.parquet"
DEF_BOUND = f"{R4}/hisss_watershed_boundaries_unsimplified.gpkg"
RATIO = 1.5
OFFICIAL = ("gagesii", "wsc_eccc")
SIMPLIFY_ABOVE_NCOORD, SIMPLIFY_TOL_M = 200000, 100.0


def key(s):
    s = str(s).strip()
    return s.lstrip("0") if s.isdigit() else s.upper()


def sha(p, blk=1 << 24):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(blk), b""):
            h.update(b)
    return h.hexdigest()


def gpkg_blob_to_geom(b):
    """GeoPackage binary -> shapely geometry (skip the GP header and its envelope)."""
    env = (b[3] >> 1) & 7
    return shapely.from_wkb(bytes(b[8 + {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[env]:]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--basinatlas", default=DEF_BASINATLAS)
    ap.add_argument("--attributes", default=DEF_ATTR, help="the deposited HydroATLAS attribute table (parquet)")
    ap.add_argument("--boundaries", default=DEF_BOUND, help="the deposited unsimplified boundary GeoPackage")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    to_alb = pyproj.Transformer.from_crs("EPSG:4326", "ESRI:102008", always_xy=True)

    def proj(g):
        return shapely.transform(g, lambda xy: np.column_stack(to_alb.transform(xy[:, 0], xy[:, 1])))

    # ---- 1. BasinATLAS attributes (no geometry) -> parquet for the R step
    info = pyogrio.read_info(a.basinatlas)
    attrs = pyogrio.read_dataframe(a.basinatlas, read_geometry=False)
    for c in ("HYBAS_ID", "NEXT_DOWN"):
        attrs[c] = attrs[c].astype("int64").astype(str)
    assert attrs.HYBAS_ID.is_unique and len(attrs) == info["features"]
    pq.write_table(pa.Table.from_pandas(attrs, preserve_index=False), os.path.join(a.out, "basinatlas_attributes.parquet"),
                   compression="snappy")
    print(f"BasinATLAS attributes: {len(attrs):,} units x {attrs.shape[1]} columns ({time.time() - t0:.0f}s)", flush=True)

    # ---- 2. level-12 polygons in ESRI:102008, drainage forest, cumulative projected area
    hb = attrs.HYBAS_ID.to_numpy(); n = len(hb)
    pos = pd.Series(np.arange(n), index=hb)
    nd = pos.reindex(attrs.NEXT_DOWN.to_numpy()).fillna(-1).astype(int).to_numpy()
    SUB = attrs.SUB_AREA.to_numpy(dtype=float)
    geoms = np.empty(n, dtype=object)
    CH = 30000
    for s in range(0, n, CH):
        df = pyogrio.read_dataframe(a.basinatlas, columns=["HYBAS_ID"], skip_features=s, max_features=CH)
        assert (df.HYBAS_ID.astype("int64").astype(str).to_numpy() == hb[s:s + len(df)]).all(), "row order differs between reads"
        g = proj(np.asarray(df.geometry.values))
        bad = ~shapely.is_valid(g)
        g[bad] = shapely.make_valid(g[bad])
        geoms[s:s + len(df)] = g
    area = shapely.area(geoms) / 1e6
    tree = shapely.STRtree(geoms)
    kids = [[] for _ in range(n)]
    for i, j in enumerate(nd):
        if j >= 0:
            kids[j].append(i)
    has_kids = np.array([len(k) > 0 for k in kids])
    tin = np.zeros(n, int); tout = np.zeros(n, int); t = 0
    for r in np.flatnonzero(nd < 0):
        st = [(r, 0)]
        while st:
            v, k = st.pop()
            if k == 0:
                tin[v] = t; t += 1
            if k < len(kids[v]):
                st.append((v, k + 1)); st.append((kids[v][k], 0))
            else:
                tout[v] = t - 1
    assert t == n, "the NEXT_DOWN graph is not a forest"
    cum = area.copy()
    for v in np.argsort(tin)[::-1]:
        if nd[v] >= 0:
            cum[nd[v]] += cum[v]
    print(f"level-12 polygons projected, forest built ({time.time() - t0:.0f}s)", flush=True)

    # ---- 3. gages: attribute table x official outlines
    A = pq.read_table(a.attributes, columns=["gage_id", "Downstream_HB_ID", "basin_area"]).to_pandas()
    A["gage_id"] = A.gage_id.astype(str); A["canon_id"] = A.gage_id.map(key)
    assert A.canon_id.is_unique
    A["c"] = pos.reindex(A.Downstream_HB_ID.astype(str).to_numpy()).to_numpy()
    assert A.c.notna().all(), "an attribute-table outlet is absent from BasinATLAS"
    A["c"] = A.c.astype(int)
    con = sqlite3.connect(f"file:{a.boundaries}?mode=ro&immutable=1", uri=True)
    L = pd.read_sql_query("select fid, canon_id, watershed_geom_source as source, geom_area_km2 from watersheds", con)
    L["canon_id"] = L.canon_id.astype(str).map(key)
    assert L.canon_id.is_unique
    A = A.merge(L, on="canon_id", how="left")
    A["official"] = A.source.isin(OFFICIAL)
    A["official_area_km2"] = A.geom_area_km2.where(A.official)
    assert (A.loc[A.official, "official_area_km2"] > 0).all()
    A["c_SUB_AREA"] = SUB[A.c]
    A["c_is_headwater"] = ~has_kids[A.c]
    A["rule_fires"] = A.official & ~A.c_is_headwater & (A.c_SUB_AREA > RATIO * A.official_area_km2)
    A["assignment_method"] = np.where(A.rule_fires, "containing_subbasin_only", "upstream")

    # ---- 4. topology check against the agency polygon
    def intersection_areas(P, cand):
        """km2 of P inside each candidate unit; complex polygons (> 20,000 vertices) are cut into 20 km tiles first
        (same areas, far faster), as in the 2026-11-19 review overlay."""
        if not len(cand):
            return np.zeros(0)
        if shapely.get_num_coordinates(P) <= 20000:
            return shapely.area(shapely.intersection(geoms[cand], P)) / 1e6
        x0, y0, x1, y1 = P.bounds; s = 20000.0
        tiles = [shapely.clip_by_rect(P, x, y, x + s, y + s) for x in np.arange(x0, x1, s) for y in np.arange(y0, y1, s)]
        tiles = np.array([t for t in tiles if not t.is_empty], dtype=object)
        a_i, t_i = shapely.STRtree(tiles).query(geoms[cand], predicate="intersects")
        ar = shapely.area(shapely.intersection(geoms[cand][a_i], tiles[t_i])) / 1e6
        return np.bincount(a_i, weights=ar, minlength=len(cand))

    rec = {k: np.full(len(A), np.nan) for k in ("P_proj_km2", "c_proj_km2", "upstream_proj_km2", "inter_c_km2", "inter_upstream_km2")}
    for k, (ix, r) in enumerate(A[A.official].iterrows()):
        blob = con.execute("select geom from watersheds where fid=?", (int(r.fid),)).fetchone()[0]
        P = proj(gpkg_blob_to_geom(blob))
        if not P.is_valid:
            P = shapely.make_valid(P)
        if P.geom_type == "GeometryCollection":
            parts = shapely.get_parts(P)
            P = shapely.union_all(parts[shapely.area(parts) > 0])
        if shapely.get_num_coordinates(P) > SIMPLIFY_ABOVE_NCOORD:
            P = shapely.simplify(P, SIMPLIFY_TOL_M, preserve_topology=True)
        cand = tree.query(P, predicate="intersects")
        ia = intersection_areas(P, cand)
        c = r.c
        inS = (tin[cand] >= tin[c]) & (tin[cand] <= tout[c])
        rec["P_proj_km2"][ix] = shapely.area(P) / 1e6
        rec["c_proj_km2"][ix] = area[c]
        rec["upstream_proj_km2"][ix] = cum[c]
        rec["inter_c_km2"][ix] = ia[cand == c].sum()
        rec["inter_upstream_km2"][ix] = ia[inS].sum()
        if (k + 1) % 2000 == 0:
            print(f"  overlay {k + 1:,}/{int(A.official.sum()):,} ({time.time() - t0:.0f}s)", flush=True)
    for kk, v in rec.items():
        A[kk] = v
    P, ac, up, ic, iu, S = (A[k].to_numpy() for k in ("P_proj_km2", "c_proj_km2", "upstream_proj_km2", "inter_c_km2",
                                                      "inter_upstream_km2", "c_SUB_AREA"))
    with np.errstate(invalid="ignore", divide="ignore"):
        exc0 = ((up - iu) - (ac - ic)) / P                       # c + upstream: area outside P, c excused
        def0 = np.maximum((P - iu) / P, 0.0)
        bad0 = (exc0 > 0.5) | ((def0 > 0.5) & ~(P < S))
        exc1 = np.where(A.rule_fires, 0.0, exc0)                  # c alone: nothing upstream
        def1 = np.where(A.rule_fires, np.maximum((P - ic) / P, 0.0), def0)
        missed_up = np.maximum(iu - ic, 0.0) / np.maximum(P - ic, 1e-9)
        exempt1 = np.where(A.rule_fires, (P < S) & (missed_up <= 0.5), P < S)
        bad1 = (exc1 > 0.5) | ((def1 > 0.5) & ~exempt1)
    off = A.official.to_numpy()
    A["exc_up_before"], A["deficit_before"], A["topology_defect_before"] = np.where(off, exc0, np.nan), np.where(off, def0, np.nan), off & bad0
    A["exc_up_after"], A["deficit_after"], A["topology_defect_after"] = np.where(off, exc1, np.nan), np.where(off, def1, np.nan), off & bad1
    A["assignment_flag"] = np.select([~A.official, A.topology_defect_after], ["no_official_polygon", "topology_defect"], "ok")

    cols = ["canon_id", "gage_id", "Downstream_HB_ID", "basin_area", "source", "official_area_km2", "c_SUB_AREA", "c_is_headwater",
            "rule_fires", "assignment_method", "assignment_flag", "P_proj_km2", "c_proj_km2", "upstream_proj_km2", "inter_c_km2",
            "inter_upstream_km2", "exc_up_before", "deficit_before", "topology_defect_before", "exc_up_after", "deficit_after",
            "topology_defect_after"]
    out = A[cols].rename(columns={"source": "watershed_geom_source"})
    out.to_csv(os.path.join(a.out, "hydroatlas_assignment_audit.csv"), index=False, float_format="%.6f")

    counts = {
        "gages": int(len(A)), "official_polygon": int(A.official.sum()),
        "no_official_polygon": int((~A.official).sum()),
        "rule_fires": int(A.rule_fires.sum()),
        "over_ratio_headwater_left_as_is": int((A.official & A.c_is_headwater & (A.c_SUB_AREA > RATIO * A.official_area_km2)).sum()),
        "topology_defects_before": int(A.topology_defect_before.sum()), "topology_defects_after": int(A.topology_defect_after.sum()),
        "fixed": int((A.topology_defect_before & ~A.topology_defect_after).sum()),
        "newly_broken": int((~A.topology_defect_before & A.topology_defect_after).sum()),
        "assignment_flag": {k: int(v) for k, v in A.assignment_flag.value_counts().items()},
    }
    prov = {"tool": "docs/benchmarks/hydroatlas_assignment_audit.py", "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "ratio": RATIO, "official_sources": list(OFFICIAL), "crs": "ESRI:102008",
            "simplify": {"above_ncoord": SIMPLIFY_ABOVE_NCOORD, "tolerance_m": SIMPLIFY_TOL_M},
            "inputs": {k: {"path": p, "bytes": os.path.getsize(p), "sha256": sha(p)}
                       for k, p in (("basinatlas", a.basinatlas), ("attributes", a.attributes), ("boundaries", a.boundaries))},
            "versions": {"python": sys.version.split()[0], "shapely": shapely.__version__, "pyproj": pyproj.__version__,
                         "pandas": pd.__version__, "pyarrow": pa.__version__},
            "counts": counts}
    with open(os.path.join(a.out, "hydroatlas_assignment_audit_provenance.json"), "w") as f:
        json.dump(prov, f, indent=2)
    print(json.dumps(counts, indent=2))
    print(f"done in {time.time() - t0:.0f}s -> {a.out}")


if __name__ == "__main__":
    main()
