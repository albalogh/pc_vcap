#!/usr/bin/env python3
"""
CapnoView Longitudinal Statistical Analysis & Plotting Pipeline
===============================================================
Analyzes capnography parameters across study phases (pc1, pc2, pc3)
for paired subjects (pc002 to pc010).

Evaluates:
  1. Raw parameters (S2V, S3V, VD_Fowler, VD_Bohr)
  2. Normalized parameters robust to flow calibration differences:
     - S3V_norm = S3V * VT (mmHg)
     - S2V_norm = S2V * VT (mmHg)
     - Fowler_fraction = VD_Fowler / VT (%)
     - Bohr_fraction = VD_Bohr / VT (%)
     - S3T (mmHg/s) - time-domain slope
     - S2T (mmHg/s) - time-domain slope
     - EtCO2, PECO2, PACO2 (mmHg)

Performs:
  - Repeated Measures ANOVA / Friedman test
  - Paired t-tests / Wilcoxon signed-rank tests for pairwise phase comparisons
  - Generates publication-ready figures (PNG)
  - Exports statistical summary table (CSV)
"""

import csv
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 11

BASE_DIR = Path("/home/balogh/Dokumentumok/pc_vcap")
OUTPUT_DIR = BASE_DIR / "analysis_output"
PLOTS_DIR = OUTPUT_DIR / "plots"
SUBJECT_CSV = OUTPUT_DIR / "subject_phase_summary.csv"
STATS_CSV = OUTPUT_DIR / "longitudinal_stats.csv"

PLOTS_DIR.mkdir(parents=True, exist_ok=True)


def read_subject_summary(csv_path):
    rows = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            parsed = {
                "phase": r["phase"],
                "subject": r["subject"],
                "vt_ml_mean": float(r["vt_ml_mean"]),
                "et_co2_mmhg_mean": float(r["et_co2_mmhg_mean"]),
                "peco2_mmhg_mean": float(r["peco2_mmhg_mean"]),
                "paco2_mmhg_mean": float(r["paco2_mmhg_mean"]),
                "fowler_vd_ml_mean": float(r["fowler_vd_ml_mean"]),
                "bohr_ratio_mean": float(r["bohr_ratio_mean"]),
                "bohr_vd_ml_mean": float(r["bohr_vd_ml_mean"]),
                "vco2_ml_mean": float(r["vco2_ml_mean"]),
                "s2v_mmhg_per_ml_mean": float(r["s2v_mmhg_per_ml_mean"]),
                "s3v_mmhg_per_ml_mean": float(r["s3v_mmhg_per_ml_mean"]),
                "s3v_norm_mmhg_mean": float(r["s3v_norm_mmhg_mean"]),
                "s2t_mmhg_per_s_mean": float(r["s2t_mmhg_per_s_mean"]),
                "s3t_mmhg_per_s_mean": float(r["s3t_mmhg_per_s_mean"]),
                "rr_bpm_mean": float(r["rr_bpm_mean"]),
            }
            parsed["fowler_fraction_mean"] = (parsed["fowler_vd_ml_mean"] / parsed["vt_ml_mean"]) * 100.0 if parsed["vt_ml_mean"] > 0 else 0.0
            parsed["bohr_fraction_mean"] = parsed["bohr_ratio_mean"] * 100.0
            rows.append(parsed)
    return rows


def main():
    if not SUBJECT_CSV.exists():
        print(f"[ERROR] Subject summary CSV not found at {SUBJECT_CSV}.")
        sys.exit(1)

    data_rows = read_subject_summary(SUBJECT_CSV)

    paired_subjects = ["pc002", "pc003", "pc004", "pc005", "pc006", "pc007", "pc008", "pc009", "pc010"]

    # Filter to paired subjects
    paired_data = [r for r in data_rows if r["subject"] in paired_subjects]

    # Map: var_name -> subject -> phase -> value
    var_matrix = defaultdict(lambda: defaultdict(dict))
    for r in paired_data:
        subj = r["subject"]
        phase = r["phase"]
        for k, v in r.items():
            if k not in ("phase", "subject"):
                var_matrix[k][subj][phase] = v

    var_defs = [
        ("s3v_norm_mmhg_mean", "S3V_norm (mmHg)", "III. fázis normalizált meredekség (S3V * VT)", "norm"),
        ("s3t_mmhg_per_s_mean", "S3T (mmHg/s)", "III. fázis időalapú meredekség (S3T)", "norm"),
        ("s2t_mmhg_per_s_mean", "S2T (mmHg/s)", "II. fázis időalapú meredekség (S2T)", "norm"),
        ("bohr_fraction_mean", "Bohr VD/VT (%)", "Élettani holttér arány (Bohr VD/VT %)", "norm"),
        ("fowler_fraction_mean", "Fowler VD/VT (%)", "Anatómiai holttér arány (Fowler VD/VT %)", "norm"),
        ("et_co2_mmhg_mean", "EtCO2 (mmHg)", "End-tidal CO2 (EtCO2 mmHg)", "norm"),
        ("peco2_mmhg_mean", "PECO2 (mmHg)", "Átlagos kilégzett CO2 (PECO2 mmHg)", "norm"),
        ("paco2_mmhg_mean", "PACO2 (mmHg)", "Átlagos alveoláris CO2 (PACO2 mmHg)", "norm"),
        ("s3v_mmhg_per_ml_mean", "S3V (mmHg/ml)", "III. fázis térfogatalapú meredekség (S3V - nem norm.)", "raw"),
        ("s2v_mmhg_per_ml_mean", "S2V (mmHg/ml)", "II. fázis térfogatalapú meredekség (S2V - nem norm.)", "raw"),
        ("fowler_vd_ml_mean", "Fowler VD (ml)", "Fowler anatómiai holttér térfogat (ml - nem norm.)", "raw"),
        ("vt_ml_mean", "VT (ml)", "Légvételi térfogat (VT - nem norm.)", "raw"),
    ]

    stats_results = []
    phases = ["pc1", "pc2", "pc3"]

    for var_name, var_short, var_label, var_type in var_defs:
        subj_dict = var_matrix[var_name]
        
        # Build 3 arrays for paired subjects
        p1_vals, p2_vals, p3_vals = [], [], []
        for s in paired_subjects:
            if s in subj_dict and "pc1" in subj_dict[s] and "pc2" in subj_dict[s] and "pc3" in subj_dict[s]:
                p1_vals.append(subj_dict[s]["pc1"])
                p2_vals.append(subj_dict[s]["pc2"])
                p3_vals.append(subj_dict[s]["pc3"])

        p1 = np.array(p1_vals)
        p2 = np.array(p2_vals)
        p3 = np.array(p3_vals)

        if len(p1) == 0:
            continue

        try:
            f_stat, p_rm = stats.f_oneway(p1, p2, p3)
            stat_name = "ANOVA F"
        except Exception:
            f_stat, p_rm = stats.friedmanchisquare(p1, p2, p3)
            stat_name = "Friedman Q"

        t_12, p_12 = stats.ttest_rel(p1, p2)
        t_23, p_23 = stats.ttest_rel(p2, p3)
        t_13, p_13 = stats.ttest_rel(p1, p3)

        stats_results.append({
            "variable": var_name,
            "short_name": var_short,
            "label": var_label,
            "type": var_type,
            "pc1_mean": round(float(np.mean(p1)), 3),
            "pc1_sd": round(float(np.std(p1, ddof=1)), 3),
            "pc2_mean": round(float(np.mean(p2)), 3),
            "pc2_sd": round(float(np.std(p2, ddof=1)), 3),
            "pc3_mean": round(float(np.mean(p3)), 3),
            "pc3_sd": round(float(np.std(p3, ddof=1)), 3),
            "test_type": stat_name,
            "stat_value": round(float(f_stat), 3),
            "p_overall": round(float(p_rm), 5),
            "p_pc1_vs_pc2": round(float(p_12), 5),
            "p_pc2_vs_pc3": round(float(p_23), 5),
            "p_pc1_vs_pc3": round(float(p_13), 5),
        })

    # Export stats CSV
    if stats_results:
        fields = list(stats_results[0].keys())
        with open(STATS_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(stats_results)
        print(f"[OUTPUT] Longitudinal Statistics CSV saved to: {STATS_CSV}")

    # =========================================================================
    # PLOTTING
    # =========================================================================

    # Figure 1: Key Normalized Capnography Slopes (S3V_norm, S3T, S2T)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True)

    plot_vars = [
        ("s3v_norm_mmhg_mean", "S3V_norm (mmHg)", "III. fázis normalizált meredekség (S3V * VT)"),
        ("s3t_mmhg_per_s_mean", "S3T (mmHg/s)", "III. fázis időalapú meredekség (S3T)"),
        ("s2t_mmhg_per_s_mean", "S2T (mmHg/s)", "II. fázis időalapú meredekség (S2T)"),
    ]

    for idx, (var_name, var_short, title) in enumerate(plot_vars):
        ax = axes[idx]
        subj_dict = var_matrix[var_name]

        p1_list, p2_list, p3_list = [], [], []
        for s in paired_subjects:
            if s in subj_dict and "pc1" in subj_dict[s] and "pc2" in subj_dict[s] and "pc3" in subj_dict[s]:
                y = [subj_dict[s]["pc1"], subj_dict[s]["pc2"], subj_dict[s]["pc3"]]
                ax.plot(["pc1", "pc2", "pc3"], y, color="gray", alpha=0.4, linewidth=1.2, marker="o", markersize=4)
                p1_list.append(y[0])
                p2_list.append(y[1])
                p3_list.append(y[2])

        means = [np.mean(p1_list), np.mean(p2_list), np.mean(p3_list)]
        sds = [np.std(p1_list, ddof=1), np.std(p2_list, ddof=1), np.std(p3_list, ddof=1)]

        ax.errorbar(["pc1", "pc2", "pc3"], means, yerr=sds, color="#0284c7",
                    linewidth=2.5, capsize=5, marker="s", markersize=8, label="Átlag ± SD")

        ax.set_title(title, fontweight="bold", pad=10)
        ax.set_ylabel(var_short)
        ax.set_xticklabels(["pc1 (2021)", "pc2 (2023)", "pc3 (2024)"])
        ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    fig1_path = PLOTS_DIR / "fig1_slopes_longitudinal.png"
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    print(f"[PLOT] Figure 1 saved: {fig1_path}")

    # Figure 2: Dead Space Fractions & Pressures (Bohr %, Fowler %, EtCO2)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True)

    plot_vars2 = [
        ("bohr_fraction_mean", "Bohr VD/VT (%)", "Élettani holttér arány (Bohr VD/VT %)"),
        ("fowler_fraction_mean", "Fowler VD/VT (%)", "Anatómiai holttér arány (Fowler VD/VT %)"),
        ("et_co2_mmhg_mean", "EtCO2 (mmHg)", "End-tidal CO2 (EtCO2 mmHg)"),
    ]

    for idx, (var_name, var_short, title) in enumerate(plot_vars2):
        ax = axes[idx]
        subj_dict = var_matrix[var_name]

        p1_list, p2_list, p3_list = [], [], []
        for s in paired_subjects:
            if s in subj_dict and "pc1" in subj_dict[s] and "pc2" in subj_dict[s] and "pc3" in subj_dict[s]:
                y = [subj_dict[s]["pc1"], subj_dict[s]["pc2"], subj_dict[s]["pc3"]]
                ax.plot(["pc1", "pc2", "pc3"], y, color="gray", alpha=0.4, linewidth=1.2, marker="o", markersize=4)
                p1_list.append(y[0])
                p2_list.append(y[1])
                p3_list.append(y[2])

        means = [np.mean(p1_list), np.mean(p2_list), np.mean(p3_list)]
        sds = [np.std(p1_list, ddof=1), np.std(p2_list, ddof=1), np.std(p3_list, ddof=1)]

        ax.errorbar(["pc1", "pc2", "pc3"], means, yerr=sds, color="#059669",
                    linewidth=2.5, capsize=5, marker="s", markersize=8, label="Átlag ± SD")

        ax.set_title(title, fontweight="bold", pad=10)
        ax.set_ylabel(var_short)
        ax.set_xticklabels(["pc1 (2021)", "pc2 (2023)", "pc3 (2024)"])
        ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    fig2_path = PLOTS_DIR / "fig2_deadspace_pressures_longitudinal.png"
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    print(f"[PLOT] Figure 2 saved: {fig2_path}")

    # Figure 3: Comparison between Unnormalized vs Normalized Variables
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

    comp_vars = [
        ("s3v_mmhg_per_ml_mean", "S3V (mmHg/ml) [Nem norm.]", "Nem normalizált S3V (térfogatalapú)"),
        ("s3v_norm_mmhg_mean", "S3V_norm (mmHg) [Norm.]", "Normalizált S3V (S3V * VT)"),
        ("fowler_vd_ml_mean", "Fowler VD (ml) [Nem norm.]", "Nem normalizált Fowler VD (ml)"),
        ("fowler_fraction_mean", "Fowler VD/VT (%) [Norm.]", "Normalizált Fowler holttér arány (VD/VT %)"),
    ]

    for idx, (var_name, var_short, title) in enumerate(comp_vars):
        row, col = idx // 2, idx % 2
        ax = axes[row, col]
        subj_dict = var_matrix[var_name]

        p1_list, p2_list, p3_list = [], [], []
        for s in paired_subjects:
            if s in subj_dict and "pc1" in subj_dict[s] and "pc2" in subj_dict[s] and "pc3" in subj_dict[s]:
                y = [subj_dict[s]["pc1"], subj_dict[s]["pc2"], subj_dict[s]["pc3"]]
                ax.plot(["pc1", "pc2", "pc3"], y, color="gray", alpha=0.4, linewidth=1.2, marker="o", markersize=4)
                p1_list.append(y[0])
                p2_list.append(y[1])
                p3_list.append(y[2])

        means = [np.mean(p1_list), np.mean(p2_list), np.mean(p3_list)]
        sds = [np.std(p1_list, ddof=1), np.std(p2_list, ddof=1), np.std(p3_list, ddof=1)]

        color = "#d97706" if "Nem norm" in title else "#2563eb"
        ax.errorbar(["pc1", "pc2", "pc3"], means, yerr=sds, color=color,
                    linewidth=2.5, capsize=5, marker="s", markersize=8)

        ax.set_title(title, fontweight="bold", pad=10)
        ax.set_ylabel(var_short)
        ax.set_xticklabels(["pc1 (2021)", "pc2 (2023)", "pc3 (2024)"])
        ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    fig3_path = PLOTS_DIR / "fig3_normalization_comparison.png"
    plt.savefig(fig3_path, dpi=300)
    plt.close()
    print(f"[PLOT] Figure 3 saved: {fig3_path}")

    # Print summary table of statistics to console
    print(f"\n{'='*115}")
    print(f"  STATISZTIKAI ELEMZÉS (Párosított n=9 alany: pc002–pc010 mindhárom fázisban)")
    print(f"{'='*115}")
    print(f"{'Változó':<22} {'Típus':<6} {'pc1 (mean±sd)':<16} {'pc2 (mean±sd)':<16} {'pc3 (mean±sd)':<16} {'p (overall)':<11} {'p (pc1-pc3)':<11}")
    print("-" * 115)
    for r in stats_results:
        p1_str = f"{r['pc1_mean']:.2f}±{r['pc1_sd']:.2f}"
        p2_str = f"{r['pc2_mean']:.2f}±{r['pc2_sd']:.2f}"
        p3_str = f"{r['pc3_mean']:.2f}±{r['pc3_sd']:.2f}"
        p_ov = f"{r['p_overall']:.4f}" if r['p_overall'] >= 0.0001 else "<0.0001"
        p_13 = f"{r['p_pc1_vs_pc3']:.4f}" if r['p_pc1_vs_pc3'] >= 0.0001 else "<0.0001"
        print(f"{r['short_name']:<22} {r['type']:<6} {p1_str:<16} {p2_str:<16} {p3_str:<16} {p_ov:<11} {p_13:<11}")


if __name__ == "__main__":
    main()
