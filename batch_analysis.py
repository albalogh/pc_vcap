#!/usr/bin/env python3
"""
CapnoView Batch Volumetric Capnography Analysis Pipeline
=========================================================
Longitudinal analysis of capnography parameters across study phases (pc1, pc2, pc3).

Extracts per-breath and per-subject parameters:
  - S2T, S2V: Phase II slope (time-based and volume-based)
  - S3T, S3V: Phase III slope (time-based and volume-based)
  - Normalized Phase III slope: S3V_norm = S3V * VT, S3T_norm = S3T / PACO2
  - Fowler anatomical dead space (VD_Fowler)
  - Bohr physiological dead space (VD_Bohr, VD/VT ratio)
  - PECO2, PACO2, EtCO2, VCO2
  - VT, RR, Ti, Te

Outputs:
  - Per-breath CSV: all breaths from all files with phase/subject/file metadata
  - Per-subject CSV: subject-level averages per phase (for longitudinal stats)
"""

import struct
import csv
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np

# Cross-version trapezoid helper (numpy 2.x support)
trapezoid = getattr(np, "trapezoid", getattr(np, "trapz", None))

# =============================================================================
# Configuration
# =============================================================================

BASE_DIR = Path("/home/balogh/Dokumentumok/pc_vcap")

# Phase definitions with calibrated gains per measurement system
PHASES = {
    "pc1": {
        "dir": BASE_DIR / "pc1",
        "dur_filter": 30.0,
        "label": "Phase 1 (2021-2022)",
        "flow_gain": 60.0,
        "co2_gain": 0.231,
        "co2_zero": 22.0,
    },
    "pc2": {
        "dir": BASE_DIR / "pc2",
        "dur_filter": 30.0,
        "label": "Phase 2 (2023)",
        "flow_gain": 18.0,
        "co2_gain": 0.231,
        "co2_zero": 22.0,
    },
    "pc3": {
        "dir": BASE_DIR / "pc3",
        "dur_filter": None,
        "label": "Phase 3 (2024)",
        "flow_gain": 19.56,
        "co2_gain": 0.231,
        "co2_zero": 22.0,
    },
}

EXCLUDED_SUBJECTS = {"teszt2", "vcal", "teszt"}

OUTPUT_DIR = BASE_DIR / "analysis_output"
BREATH_CSV = OUTPUT_DIR / "all_breaths.csv"
SUBJECT_CSV = OUTPUT_DIR / "subject_phase_summary.csv"


# =============================================================================
# .inp Decoder
# =============================================================================

def read_pascal_u16_string(blob: bytes, pos: int) -> tuple:
    length = struct.unpack_from("<H", blob, pos)[0]
    pos += 2
    value = blob[pos:pos + length].decode("ascii", errors="replace")
    return value, pos + length


def decode_inp(path: Path) -> dict:
    blob = path.read_bytes()
    if len(blob) < 44:
        raise ValueError("File too short")

    header = {
        "sample_rate_hz": float(struct.unpack_from("<f", blob, 0)[0]),
        "sample_count": int(struct.unpack_from("<H", blob, 4)[0]),
        "channel_count": int(struct.unpack_from("<H", blob, 6)[0]),
        "duration_s": float(struct.unpack_from("<f", blob, 8)[0]),
    }

    pos = 42
    header["recorded_at"], pos = read_pascal_u16_string(blob, pos)
    header["subject_id"], pos = read_pascal_u16_string(blob, pos)

    if pos + 2 > len(blob):
        raise ValueError("Truncated metadata")
    pos += 2

    header["signal_configuration_path"], pos = read_pascal_u16_string(blob, pos)

    payload = np.frombuffer(blob, dtype="<i2", offset=pos)
    expected = header["sample_count"] * header["channel_count"]
    if payload.size != expected:
        raise ValueError(f"Payload mismatch: {payload.size} vs {expected}")

    signals = payload.reshape(header["sample_count"], header["channel_count"])
    return header, signals


# =============================================================================
# Signal Processing & Capnography Parameter Extraction
# =============================================================================

def process_file(path: Path, phase: str, phase_cfg: dict) -> list:
    header, signals = decode_inp(path)
    sr = header["sample_rate_hz"]
    n = header["sample_count"]
    dt = 1.0 / sr

    ch1_raw = signals[:, 0].astype(np.float64)
    ch4_raw = signals[:, 3].astype(np.float64) if header["channel_count"] >= 4 else None

    flow_gain = phase_cfg["flow_gain"]
    co2_gain = phase_cfg["co2_gain"]
    co2_zero = phase_cfg["co2_zero"]

    # Flow calibration
    flow_base = float(np.median(ch1_raw))
    flow_mls = -(ch1_raw - flow_base) * flow_gain

    # CO2 calibration
    ch4_arr = None
    if ch4_raw is not None:
        ch4_arr = (ch4_raw - co2_zero) * co2_gain

    time_s = np.arange(n) * dt

    # Robust Breath Detection
    breath_starts = [0]
    use_co2 = False

    if ch4_arr is not None:
        min_co2 = float(np.min(ch4_arr))
        max_co2 = float(np.max(ch4_arr))
        if (max_co2 - min_co2) > 5.0:
            use_co2 = True
            thresh_high = min_co2 + (max_co2 - min_co2) * 0.35
            thresh_low = min_co2 + (max_co2 - min_co2) * 0.20
            in_expiration = False
            for i in range(1, n):
                if ch4_arr[i] > thresh_high:
                    in_expiration = True
                elif in_expiration and ch4_arr[i] < thresh_low:
                    if (time_s[i] - time_s[breath_starts[-1]]) > 0.8:
                        breath_starts.append(i)
                        in_expiration = False

    if not use_co2 or len(breath_starts) < 3:
        breath_starts = [0]
        in_neg = False
        for i in range(1, n):
            if flow_mls[i] < -20.0:
                in_neg = True
            elif in_neg and flow_mls[i] > 20.0:
                if (time_s[i] - time_s[breath_starts[-1]]) > 0.8:
                    breath_starts.append(i)
                    in_neg = False

    if len(breath_starts) < 2:
        return []

    if breath_starts[-1] != n - 1:
        breath_starts.append(n - 1)

    results = []

    for k in range(len(breath_starts) - 1):
        i0 = breath_starts[k]
        i1 = breath_starts[k + 1]

        if (i1 - i0) < 10:
            continue

        cycle_flow = flow_mls[i0:i1 + 1]
        cycle_co2 = ch4_arr[i0:i1 + 1] if ch4_arr is not None else None
        cycle_time = time_s[i0:i1 + 1]
        cycle_dt = dt

        # Volume via cumulative integration with linear drift correction
        cycle_vol_raw = np.cumsum(cycle_flow) * cycle_dt
        drift = cycle_vol_raw[-1]
        correction = np.linspace(0, drift, len(cycle_vol_raw))
        cycle_vol = cycle_vol_raw - correction

        peak_rel = int(np.argmax(cycle_vol))
        vt = float(cycle_vol[peak_rel])

        if vt < 50:
            continue

        ti = float(cycle_time[peak_rel] - cycle_time[0])
        te = float(cycle_time[-1] - cycle_time[peak_rel])
        tot_time = ti + te
        rr = 60.0 / tot_time if tot_time > 0 else 0.0

        et_co2 = 0.0
        peco2 = 0.0
        paco2 = 0.0
        fowler_vd = 0.0
        bohr_ratio = 0.0
        bohr_vd = 0.0
        vco2_ml = 0.0
        s2v = 0.0
        s3v = 0.0
        s2t = 0.0
        s3t = 0.0

        if cycle_co2 is not None and (i1 - (i0 + peak_rel)) > 5:
            co2_exp = cycle_co2[peak_rel:]
            v_exh = vt - cycle_vol[peak_rel:]
            et_co2 = float(np.max(co2_exp))

            sort_idx = np.argsort(v_exh)
            v_sorted = v_exh[sort_idx]
            co2_sorted = co2_exp[sort_idx]
            _, uniq_idx = np.unique(v_sorted, return_index=True)
            v_uniq = v_sorted[uniq_idx]
            co2_uniq = co2_sorted[uniq_idx]

            n_grid = max(int(vt), 50)
            v_grid = np.linspace(0, vt, n_grid)
            co2_grid = np.interp(v_grid, v_uniq, co2_uniq)

            # PECO2 (volume-weighted mean expired CO2)
            peco2 = float(trapezoid(co2_grid, v_grid) / vt) if vt > 0 else 0.0

            # Phase II inflection point
            dco2_dv = np.gradient(co2_grid, v_grid)
            kernel_size = max(5, n_grid // 30)
            if kernel_size % 2 == 0:
                kernel_size += 1
            dco2_smooth = np.convolve(dco2_dv, np.ones(kernel_size) / kernel_size, mode='same')

            i_search_min = max(1, int(0.10 * n_grid))
            i_search_max = min(n_grid - 1, int(0.65 * n_grid))

            if i_search_max > i_search_min + 5:
                idx_inflect = i_search_min + int(np.argmax(dco2_smooth[i_search_min:i_search_max]))
                fowler_vd = float(v_grid[idx_inflect])
            else:
                idx_inflect = int(0.30 * n_grid)
                fowler_vd = 0.3 * vt

            # S2V (Phase II volume slope)
            s2_start = max(0, idx_inflect - int(0.10 * n_grid))
            s2_end = min(n_grid - 1, idx_inflect + int(0.10 * n_grid))
            if s2_end > s2_start + 3:
                s2v_coeffs = np.polyfit(v_grid[s2_start:s2_end + 1], co2_grid[s2_start:s2_end + 1], 1)
                s2v = float(s2v_coeffs[0])

            # PACO2 and S3V (Phase III volume slope)
            i_p3_start = min(idx_inflect + int(0.15 * n_grid), int(0.65 * n_grid))
            i_p3_end = max(i_p3_start + 5, n_grid - int(0.05 * n_grid))

            if i_p3_end > i_p3_start + 3 and i_p3_start < n_grid:
                paco2 = float(np.mean(co2_grid[i_p3_start:i_p3_end]))
                s3v_coeffs = np.polyfit(v_grid[i_p3_start:i_p3_end], co2_grid[i_p3_start:i_p3_end], 1)
                s3v = float(s3v_coeffs[0])
            else:
                paco2 = float(np.max(co2_grid))
                s3v = 0.0

            # Time-domain slopes S2T, S3T
            exp_time = cycle_time[peak_rel:] - cycle_time[peak_rel]
            dco2_dt_raw = np.gradient(co2_exp, exp_time)
            dco2_dt_smooth = np.convolve(dco2_dt_raw, np.ones(kernel_size) / kernel_size, mode='same')

            n_exp = len(exp_time)
            t_search_min = max(1, int(0.10 * n_exp))
            t_search_max = min(n_exp - 1, int(0.65 * n_exp))

            if t_search_max > t_search_min + 5:
                t_inflect = t_search_min + int(np.argmax(dco2_dt_smooth[t_search_min:t_search_max]))

                ts2_start = max(0, t_inflect - int(0.10 * n_exp))
                ts2_end = min(n_exp - 1, t_inflect + int(0.10 * n_exp))
                if ts2_end > ts2_start + 3:
                    s2t_coeffs = np.polyfit(exp_time[ts2_start:ts2_end + 1], co2_exp[ts2_start:ts2_end + 1], 1)
                    s2t = float(s2t_coeffs[0])

                tp3_start = min(t_inflect + int(0.15 * n_exp), int(0.65 * n_exp))
                tp3_end = max(tp3_start + 5, n_exp - int(0.05 * n_exp))
                if tp3_end > tp3_start + 3:
                    s3t_coeffs = np.polyfit(exp_time[tp3_start:tp3_end], co2_exp[tp3_start:tp3_end], 1)
                    s3t = float(s3t_coeffs[0])

            # Bohr dead space & VCO2
            bohr_ratio = (paco2 - peco2) / paco2 if paco2 > 0 else 0.0
            bohr_vd = bohr_ratio * vt
            vco2_ml = (peco2 / 760.0) * vt

        results.append({
            "phase": phase,
            "subject": header["subject_id"].lower(),
            "file": path.name,
            "date": header["recorded_at"].strip(),
            "breath_num": k + 1,
            "vt_ml": round(vt, 1),
            "et_co2_mmhg": round(et_co2, 1),
            "peco2_mmhg": round(peco2, 2),
            "paco2_mmhg": round(paco2, 2),
            "fowler_vd_ml": round(fowler_vd, 1),
            "bohr_ratio": round(bohr_ratio, 4),
            "bohr_vd_ml": round(bohr_vd, 1),
            "vco2_ml": round(vco2_ml, 2),
            "s2v_mmhg_per_ml": round(s2v, 4),
            "s3v_mmhg_per_ml": round(s3v, 4),
            "s3v_norm_mmhg": round(s3v * vt, 2),
            "s2t_mmhg_per_s": round(s2t, 2),
            "s3t_mmhg_per_s": round(s3t, 2),
            "ti_s": round(ti, 3),
            "te_s": round(te, 3),
            "rr_bpm": round(rr, 1),
        })

    return results


# =============================================================================
# Main Pipeline
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_breaths = []
    file_count = 0
    skipped_count = 0

    for phase_name, phase_cfg in PHASES.items():
        phase_dir = phase_cfg["dir"]
        dur_filter = phase_cfg["dur_filter"]

        if not phase_dir.exists():
            print(f"[WARNING] Phase directory not found: {phase_dir}")
            continue

        inp_files = sorted(phase_dir.glob("*.inp"), key=lambda p: p.name.lower())
        print(f"\n{'='*60}")
        print(f"  {phase_cfg['label']} ({phase_name}): {len(inp_files)} .inp files")
        print(f"{'='*60}")

        for fpath in inp_files:
            try:
                header, _ = decode_inp(fpath)
            except Exception as e:
                print(f"  [SKIP] {fpath.name}: decode error: {e}")
                skipped_count += 1
                continue

            subject = header["subject_id"].lower()

            if subject in EXCLUDED_SUBJECTS:
                continue

            if dur_filter is not None and abs(header["duration_s"] - dur_filter) > 1.0:
                continue

            try:
                breaths = process_file(fpath, phase_name, phase_cfg)
                if breaths:
                    file_count += 1
                    all_breaths.extend(breaths)
                    n_breaths = len(breaths)
                    avg_vt = np.mean([b["vt_ml"] for b in breaths])
                    avg_etco2 = np.mean([b["et_co2_mmhg"] for b in breaths])
                    print(f"  [OK] {fpath.name}: {n_breaths} breaths, avg VT={avg_vt:.0f}ml, avg EtCO2={avg_etco2:.1f}mmHg")
                else:
                    print(f"  [NO BREATHS] {fpath.name}")
            except Exception as e:
                print(f"  [ERROR] {fpath.name}: {e}")
                skipped_count += 1

    if not all_breaths:
        print("\n[ERROR] No breaths extracted. Check calibration and data.")
        return

    # Write per-breath CSV
    breath_fields = list(all_breaths[0].keys())
    with open(BREATH_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=breath_fields)
        writer.writeheader()
        writer.writerows(all_breaths)
    print(f"\n[OUTPUT] Per-breath CSV: {BREATH_CSV} ({len(all_breaths)} rows)")

    # Compute per-subject per-phase averages
    subject_phase = defaultdict(list)
    for b in all_breaths:
        key = (b["phase"], b["subject"])
        subject_phase[key].append(b)

    summary_rows = []
    numeric_keys = [
        "vt_ml", "et_co2_mmhg", "peco2_mmhg", "paco2_mmhg",
        "fowler_vd_ml", "bohr_ratio", "bohr_vd_ml", "vco2_ml",
        "s2v_mmhg_per_ml", "s3v_mmhg_per_ml", "s3v_norm_mmhg",
        "s2t_mmhg_per_s", "s3t_mmhg_per_s",
        "ti_s", "te_s", "rr_bpm"
    ]

    for (phase, subject), breaths in sorted(subject_phase.items()):
        row = {
            "phase": phase,
            "subject": subject,
            "n_files": len(set(b["file"] for b in breaths)),
            "n_breaths": len(breaths),
            "date_first": min(b["date"] for b in breaths),
            "date_last": max(b["date"] for b in breaths),
        }
        for key in numeric_keys:
            values = [b[key] for b in breaths]
            row[f"{key}_mean"] = round(float(np.mean(values)), 3)
            row[f"{key}_sd"] = round(float(np.std(values, ddof=1)), 3) if len(values) > 1 else 0.0
        summary_rows.append(row)

    summary_fields = list(summary_rows[0].keys())
    with open(SUBJECT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"[OUTPUT] Subject-phase summary CSV: {SUBJECT_CSV} ({len(summary_rows)} rows)")

    # Print summary table
    print(f"\n{'='*105}")
    print(f"  ÖSSZEFOGLALÓ — {len(summary_rows)} alany×fázis kombináció")
    print(f"{'='*105}")
    print(f"{'Phase':<6} {'Subject':<8} {'N_br':>5} {'VT':>6} {'EtCO2':>6} {'S2T':>7} {'S2V':>8} {'S3T':>7} {'S3V':>8} {'FowlVD':>7} {'Bohr%':>6}")
    print("-" * 105)
    for r in summary_rows:
        print(f"{r['phase']:<6} {r['subject']:<8} {r['n_breaths']:>5} "
              f"{r['vt_ml_mean']:>6.0f} {r['et_co2_mmhg_mean']:>6.1f} "
              f"{r['s2t_mmhg_per_s_mean']:>7.1f} {r['s2v_mmhg_per_ml_mean']:>8.4f} "
              f"{r['s3t_mmhg_per_s_mean']:>7.2f} {r['s3v_mmhg_per_ml_mean']:>8.4f} "
              f"{r['fowler_vd_ml_mean']:>7.1f} {r['bohr_ratio_mean']*100:>5.1f}%")

    print(f"\nProcessed {file_count} files, skipped {skipped_count}, extracted {len(all_breaths)} breaths total.")


if __name__ == "__main__":
    main()
