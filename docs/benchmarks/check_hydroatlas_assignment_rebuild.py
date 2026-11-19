#!/usr/bin/env python
"""GATE: the rebuilt HydroATLAS attribute table against the deposited one (containing-subbasin rule, 2026-11-19).

Checks (exit 1 on any failure):
  1. shape and schema: the same 8,014 gages in the same order; the deposited 211 columns with the same names, order
     and types, plus `assignment_method` and `assignment_flag` after `outlet_COAST`;
  2. rows with assignment_method = upstream: every one of the 211 columns is identical to the deposited parquet
     (NaN = NaN), and the CSV cells are character-identical to the deposited CSV;
  3. rows with assignment_method = containing_subbasin_only: every attribute is re-derived here, independently of the
     R code, from the containing subbasin's own BasinATLAS values (_s twin, -9999 NoData -> NA; NA where no
     subbasin-level value exists; class argmax over the _s fractions; single-member area diagnostics);
  4. the rule decision and the flags agree with the audit table, and the expected counts hold;
  5. the CSV and the parquet of the rebuild agree (floats to 1e-8), and the dictionary keeps every deposited row's
     column / theme / unit, adding the two new columns.
Usage:
  python docs/benchmarks/check_hydroatlas_assignment_rebuild.py --new <dir> --audit-dir <dir> [--deposited <dir>]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HOME = os.path.expanduser("~")
DEF_DEP = f"{HOME}/Downloads/Signatures_v3_28oct2026/resource4_geometry_attributes"
DEF_R3 = f"{HOME}/Downloads/Signatures_v3_28oct2026/resource3_input_data"
ATTR = "hisss_watershed_hydroatlas_attributes"
DICT = "hisss_hydroatlas_attributes_dictionary.csv"
NEW_COLS = ["assignment_method", "assignment_flag"]
NO_TWIN = {"lkv_mc_usu", "rev_mc_usu", "dis_m3_pyr", "dis_m3_pmn", "dis_m3_pmx", "dor_pc_pva"}
FRACS = {"glc_cl_smj": [f"glc_pc_s{i:02d}" for i in range(1, 23)], "pnv_cl_smj": [f"pnv_pc_s{i:02d}" for i in range(1, 16)]}


def key(s):
    s = str(s).strip()
    return s.lstrip("0") if s.isdigit() else s.upper()


def same(a, b):
    a = pd.Series(a).reset_index(drop=True); b = pd.Series(b).reset_index(drop=True)
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        x, y = a.to_numpy(dtype=float), b.to_numpy(dtype=float)
        return (x == y) | (np.isnan(x) & np.isnan(y))
    return (a.astype(str) == b.astype(str)).to_numpy()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", required=True)
    ap.add_argument("--audit-dir", required=True)
    ap.add_argument("--deposited", default=DEF_DEP)
    a = ap.parse_args()
    fails = []

    def check(ok, msg):
        print(("PASS " if ok else "FAIL ") + msg)
        if not ok:
            fails.append(msg)

    dep_t = pq.read_table(os.path.join(a.deposited, ATTR + ".parquet"))
    new_t = pq.read_table(os.path.join(a.new, ATTR + ".parquet"))
    dep, new = dep_t.to_pandas(), new_t.to_pandas()
    dcols = list(dep.columns)
    i_coast = dcols.index("outlet_COAST") + 1
    check(list(new.columns) == dcols[:i_coast] + NEW_COLS + dcols[i_coast:], f"columns: {len(dcols)} deposited + {NEW_COLS} after outlet_COAST")
    dtyp = {f.name: str(f.type) for f in dep_t.schema}; ntyp = {f.name: str(f.type) for f in new_t.schema}
    bad_t = [c for c in dcols if dtyp[c] != ntyp.get(c)]
    check(not bad_t, f"types of the {len(dcols)} deposited columns unchanged {bad_t[:5]}")
    check(len(new) == len(dep) and (new.gage_id.astype(str).to_numpy() == dep.gage_id.astype(str).to_numpy()).all(),
          f"same {len(dep):,} gages in the same order")

    up = (new.assignment_method == "upstream").to_numpy(); loc = (new.assignment_method == "containing_subbasin_only").to_numpy()
    check((up | loc).all(), "assignment_method is upstream or containing_subbasin_only on every row")
    diff_up = {c: int((~same(new.loc[up, c], dep.loc[up, c])).sum()) for c in dcols}
    nz = {c: v for c, v in diff_up.items() if v}
    check(not nz, f"{int(up.sum()):,} upstream rows: all {len(dcols)} columns identical to the deposited parquet {dict(list(nz.items())[:5])}")
    id_cols = ["gage_id", "latitude", "longitude", "basin_area", "gage_type", "Downstream_HB_ID"]
    check(all(same(new.loc[loc, c], dep.loc[loc, c]).all() for c in id_cols), f"{int(loc.sum())} rule rows: identifiers and basin_area unchanged")

    # ---- 3. independent re-derivation of the rule rows from BasinATLAS
    ba = pd.read_parquet(os.path.join(a.audit_dir, "basinatlas_attributes.parquet")).set_index("HYBAS_ID")
    L = new.loc[loc].copy()
    c = ba.loc[L.Downstream_HB_ID.astype(str).to_numpy()]
    c.index = L.index
    exp = {"n_member_basins": np.ones(len(L)), "watershed_area_sum_SUB_AREA": c.SUB_AREA.to_numpy(float),
           "outlet_UP_AREA": c.UP_AREA.to_numpy(float), "outlet_ORDER": c.ORDER_.to_numpy(float),
           "outlet_ENDO": c.ENDO.to_numpy(float), "outlet_COAST": c.COAST.to_numpy(float)}
    ba_area = L.basin_area.to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        exp["watershed_area_rel_diff"] = np.where(np.isfinite(ba_area) & (ba_area > 0), np.abs(c.SUB_AREA.to_numpy(float) - ba_area) / ba_area, np.nan)
    attr_cols = dcols[dcols.index("dis_m3_pyr"):]
    for col in attr_cols:
        if col.endswith("_frac") or col in FRACS:
            base = col[:-5] if col.endswith("_frac") else col
            fr = c[FRACS[base]].to_numpy(float).copy(); fr[fr == -9999] = np.nan
            mx = np.nanmax(np.where(np.isnan(fr), -np.inf, fr), axis=1)
            code = np.array([int(FRACS[base][i][-2:]) for i in np.argmax(np.where(np.isnan(fr), -np.inf, fr), axis=1)], float)
            ok = np.isfinite(mx) & (mx > 0)
            exp[col] = np.where(ok, mx, np.nan) if col.endswith("_frac") else np.where(ok, code, np.nan)
            continue
        ext = col.split("_")[2][0]
        if col in NO_TWIN or ext == "p":
            exp[col] = np.full(len(L), np.nan)
            continue
        src = col[:7] + "s" + col[8:] if ext == "u" else col
        v = c[src].to_numpy(float)
        unit = col.split("_")[1]
        if unit == "cl" or col == "gad_id_smj":
            v = np.where(v < 0, np.nan, v)
        elif col not in ("ele_mt_smn", "ele_mt_smx"):
            v = np.where(v == -9999, np.nan, v)
        exp[col] = v
    bad = {}
    for col, e in exp.items():
        got = L[col].to_numpy(float)
        m = ~((got == e) | (np.isnan(got) & np.isnan(e)))
        if col == "watershed_area_rel_diff":
            m = ~((np.abs(got - e) <= 1e-12) | (np.isnan(got) & np.isnan(e)))
        if m.any():
            bad[col] = int(m.sum())
    check(not bad, f"{len(L)} rule rows: all {len(exp)} derived columns match an independent BasinATLAS re-derivation {dict(list(bad.items())[:6])}")
    changed_cols = [c_ for c_ in dcols if (~same(L[c_], dep.loc[loc, c_])).any()]
    print(f"INFO rule rows: {len(changed_cols)} of {len(dcols)} columns differ from the deposited values")

    # ---- 4. audit agreement and expected counts
    au = pd.read_csv(os.path.join(a.audit_dir, "hydroatlas_assignment_audit.csv"), dtype={"canon_id": str, "gage_id": str, "Downstream_HB_ID": str})
    au = au.set_index("canon_id").loc[new.gage_id.astype(str).map(key)]
    check((au.rule_fires.astype(str).str.upper().eq("TRUE").to_numpy() == loc).all(), "rule decision identical to the audit table")
    check((au.assignment_flag.to_numpy() == new.assignment_flag.to_numpy()).all(), "assignment_flag identical to the audit table")
    sel = {w: set(map(key, open(os.path.join(DEF_R3, f"hisss_gage_selection_wy{w}_2025.txt")).read().split())) for w in ("1993", "1980")}
    canon = new.gage_id.astype(str).map(key)
    n_loc = {w: int((loc & canon.isin(s).to_numpy()).sum()) for w, s in sel.items()}
    fl = new.assignment_flag.value_counts().to_dict()
    print(f"INFO containing_subbasin_only {int(loc.sum())} (WY1993-2025 {n_loc['1993']}, WY1980-2025 {n_loc['1980']}); flags {fl}")
    check(int(loc.sum()) == 511 and n_loc == {"1993": 424, "1980": 352}, "511 rule rows (424 / 352 in the two products)")
    check(fl.get("topology_defect") == 321 and fl.get("no_official_polygon") == 31, "flags: 321 topology_defect, 31 no_official_polygon")

    # ---- 5. CSV vs parquet (rebuild) and deposited CSV cells on upstream rows; dictionary
    ncsv = pd.read_csv(os.path.join(a.new, ATTR + ".csv"), dtype=str, keep_default_na=False)
    dcsv = pd.read_csv(os.path.join(a.deposited, ATTR + ".csv"), dtype=str, keep_default_na=False)
    check(list(ncsv.columns) == list(new.columns), "rebuild CSV columns = parquet columns")
    cell = {c_: int((ncsv.loc[up, c_].to_numpy() != dcsv.loc[up, c_].to_numpy()).sum()) for c_ in dcols}
    nzc = {k: v for k, v in cell.items() if v}
    check(not nzc, f"{int(up.sum()):,} upstream rows: CSV cells character-identical to the deposited CSV {dict(list(nzc.items())[:5])}")
    mx = 0.0
    for c_ in new.columns:
        if pd.api.types.is_numeric_dtype(new[c_]):
            x = pd.to_numeric(ncsv[c_].replace("", np.nan), errors="coerce").to_numpy(float); y = new[c_].to_numpy(float)
            if not (np.isnan(x) == np.isnan(y)).all():
                fails.append(f"CSV/parquet NA pattern differs in {c_}"); print(f"FAIL CSV/parquet NA pattern differs in {c_}")
            elif (~np.isnan(x)).any():
                mx = max(mx, float(np.nanmax(np.abs(x - y))))
        else:
            if not (ncsv[c_].to_numpy() == new[c_].fillna("").astype(str).to_numpy()).all():
                fails.append(f"CSV/parquet strings differ in {c_}"); print(f"FAIL CSV/parquet strings differ in {c_}")
    check(mx < 1e-8, f"rebuild CSV and parquet agree (max abs float difference {mx:.1e})")
    dd = pd.read_csv(os.path.join(a.deposited, DICT), dtype=str, keep_default_na=False)
    nd = pd.read_csv(os.path.join(a.new, DICT), dtype=str, keep_default_na=False)
    check(nd.column.tolist() == list(new.columns), "dictionary has one row per column, in table order")
    nd2 = nd.set_index("column").loc[dd.column]
    check(all((nd2[k].to_numpy() == dd[k].to_numpy()).all() for k in ("theme", "unit")), "dictionary theme / unit unchanged for every deposited column")
    ch = {k: int((nd2[k].to_numpy() != dd[k].to_numpy()).sum()) for k in ("aggregation", "description")}
    print(f"INFO dictionary rows with revised text: {ch}")

    print(f"\n{'ALL CHECKS PASS' if not fails else str(len(fails)) + ' CHECK(S) FAILED'}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
