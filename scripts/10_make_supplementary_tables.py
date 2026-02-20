#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

def must_read_csv(p: Path) -> pd.DataFrame:
    if not p.exists():
        raise SystemExit(f"NOT FOUND: {p}")
    return pd.read_csv(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reports-dir", default="reports")
    ap.add_argument("--outdir", default="submission_supplement")
    args = ap.parse_args()

    reports = Path(args.reports_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------- Supplementary Table S1 ----------
    # Combine q=0.85, 0.90, 0.95 correlation summaries
    q_files = [
        (0.85, reports / "b_th2_th17_correlations_q085.csv"),
        (0.90, reports / "b_th2_th17_correlations.csv"),
        (0.95, reports / "b_th2_th17_correlations_q095.csv"),
    ]
    dfs = []
    for q, fp in q_files:
        df = must_read_csv(fp).copy()
        df["quantile_q"] = q
        dfs.append(df)

    s1 = pd.concat(dfs, ignore_index=True)

    # Normalize/standardize expected columns
    # Keep all columns, but enforce a canonical ordering for key outputs.
    key_cols = [
        "section", "quantile_q", "n_spots", "immune_only", "immune_q",
        "spearman_B_vs_Th2", "spearman_B_vs_Th17",
        "pearson_B_vs_Th2", "pearson_B_vs_Th17",
        "spearman_adj_B_vs_Th2", "spearman_adj_B_vs_Th17",
        "perm_spearman_adj_B_vs_Th2", "perm_p_adj_B_vs_Th2",
        "perm_spearman_adj_B_vs_Th17", "perm_p_adj_B_vs_Th17",
    ]
    cols = [c for c in key_cols if c in s1.columns] + [c for c in s1.columns if c not in key_cols]
    s1 = s1[cols].sort_values(["quantile_q", "section"])

    s1_path = outdir / "Supplementary_Table_S1_correlations_quantile_sensitivity.csv"
    s1.to_csv(s1_path, index=False)
    print("[wrote]", s1_path)

    # Optional XLSX
    try:
        xlsx_path = outdir / "Supplementary_Tables.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as xw:
            s1.to_excel(xw, sheet_name="TableS1_Correlations", index=False)
        print("[wrote]", xlsx_path)
    except Exception as e:
        print("[warn] Could not write XLSX:", e)

    # ---------- Supplementary Table S2 ----------
    # LGALS9 gradient analysis (use whatever is present)
    lg_files = sorted(reports.glob("lgals9_distance_effects*.csv"))
    if not lg_files:
        raise SystemExit("No LGALS9 gradient result files found (expected reports/lgals9_distance_effects*.csv).")

    lg_dfs = []
    for fp in lg_files:
        df = pd.read_csv(fp).copy()
        df["source_file"] = fp.name
        lg_dfs.append(df)

    s2 = pd.concat(lg_dfs, ignore_index=True)

    # Sort by section then file
    if "section" in s2.columns:
        s2 = s2.sort_values(["section", "source_file"])
    else:
        s2 = s2.sort_values(["source_file"])

    s2_path = outdir / "Supplementary_Table_S2_lgals9_distance_gradient.csv"
    s2.to_csv(s2_path, index=False)
    print("[wrote]", s2_path)

    # Add to XLSX if possible
    try:
        xlsx_path = outdir / "Supplementary_Tables.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as xw:
            s2.to_excel(xw, sheet_name="TableS2_LGALS9", index=False)
        print("[updated]", xlsx_path)
    except Exception as e:
        print("[warn] Could not append Table S2 to XLSX:", e)

if __name__ == "__main__":
    main()

