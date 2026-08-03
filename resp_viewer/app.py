#!/usr/bin/env python3
"""
CapnoView Backend Server
A lightweight HTTP server and API for decoding .inp files and .csv recordings,
performing respiratory flow integration with end-expiratory detrending,
and serving an interactive multi-track medical signal viewer with advanced
volumetric capnography analysis (Fowler anatomical dead space, PECO2, PACO2, Bohr physiological dead space, VCO2).
Supports post-hoc calibration for mainstream capnograph and pneumotachograph airflow systems.
"""

import os
import sys
import json
import struct
import argparse
from pathlib import Path
from http.server import HTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

import numpy as np

WORKSPACE_DIR = Path(__file__).resolve().parent.parent


def read_pascal_u16_string(blob: bytes, pos: int) -> tuple[str, int]:
    if pos + 2 > len(blob):
        raise ValueError("Truncated string length.")
    length = struct.unpack_from("<H", blob, pos)[0]
    pos += 2
    if pos + length > len(blob):
        raise ValueError("Truncated string content.")
    value = blob[pos:pos + length].decode("ascii", errors="replace")
    return value, pos + length


def decode_inp_file(path: Path) -> dict:
    blob = path.read_bytes()
    if len(blob) < 44:
        raise ValueError("File is too short.")

    header = {
        "sample_rate_hz": float(struct.unpack_from("<f", blob, 0)[0]),
        "sample_count": int(struct.unpack_from("<H", blob, 4)[0]),
        "channel_count": int(struct.unpack_from("<H", blob, 6)[0]),
        "duration_s": float(struct.unpack_from("<f", blob, 8)[0]),
        "marker": int(struct.unpack_from("<i", blob, 36)[0]),
        "flags": int(struct.unpack_from("<H", blob, 40)[0]),
    }

    pos = 42
    header["recorded_at"], pos = read_pascal_u16_string(blob, pos)
    header["subject_id"], pos = read_pascal_u16_string(blob, pos)

    if pos + 2 > len(blob):
        raise ValueError("Truncated metadata.")
    header["unknown_u16"] = int(struct.unpack_from("<H", blob, pos)[0])
    pos += 2

    header["signal_configuration_path"], pos = read_pascal_u16_string(blob, pos)
    header["data_offset"] = pos

    payload = np.frombuffer(blob, dtype="<i2", offset=pos)
    expected = header["sample_count"] * header["channel_count"]
    if payload.size != expected:
        raise ValueError(
            f"Payload has {payload.size} int16 values; expected {expected}."
        )

    signals = payload.reshape(header["sample_count"], header["channel_count"])
    time_s = (np.arange(header["sample_count"]) / header["sample_rate_hz"]).tolist()
    
    channels = {}
    for idx in range(header["channel_count"]):
        channels[f"raw_ch{idx + 1}"] = signals[:, idx].tolist()

    return {
        "filename": path.name,
        "file_type": "inp",
        "header": header,
        "time_s": time_s,
        "channels": channels,
    }


def decode_csv_file(path: Path) -> dict:
    import csv
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        headers = next(reader)
        data = []
        for row in reader:
            if row and len(row) == len(headers):
                try:
                    data.append([float(val) for val in row])
                except ValueError:
                    continue

    arr = np.array(data)
    time_s = arr[:, 0].tolist() if arr.shape[1] > 0 else []
    sample_count = len(time_s)
    sample_rate_hz = 1.0 / (time_s[1] - time_s[0]) if sample_count > 1 and time_s[1] > time_s[0] else 256.0
    duration_s = time_s[-1] if sample_count > 0 else 0.0

    channels = {}
    for col_idx in range(1, arr.shape[1]):
        col_name = headers[col_idx]
        channels[col_name] = arr[:, col_idx].tolist()
        if "ch1" in col_name.lower() or "respiratory" in col_name.lower():
            channels["raw_ch1"] = arr[:, col_idx].tolist()
        if "ch4" in col_name.lower() or "capnogram" in col_name.lower() or "co2" in col_name.lower():
            channels["raw_ch4"] = arr[:, col_idx].tolist()

    header = {
        "sample_rate_hz": float(sample_rate_hz),
        "sample_count": sample_count,
        "channel_count": arr.shape[1] - 1,
        "duration_s": float(duration_s),
        "subject_id": path.stem,
        "recorded_at": "CSV export",
    }

    return {
        "filename": path.name,
        "file_type": "csv",
        "header": header,
        "time_s": time_s,
        "channels": channels,
    }


def compute_volume_and_breaths(time_s, ch1, ch4=None, mode="linear_detrend", polarity=1.0,
                               flow_offset=None, flow_gain=19.56, co2_zero=22.0, co2_gain=0.231):
    """
    Server-side helper to compute calibrated detrended volume signals, breath statistics,
    and volumetric capnography parameters (Fowler VD, PECO2, PACO2, Bohr VD, VCO2).
    """
    time_arr = np.array(time_s)
    ch1_arr = np.array(ch1)
    n = len(time_arr)
    if n == 0:
        return {}

    dt = time_arr[1] - time_arr[0] if n > 1 else 1.0 / 256.0

    # 1. Flow baseline & calibration (ml/s)
    base = flow_offset if flow_offset is not None else np.median(ch1_arr)
    flow_arr = -(ch1_arr - base) * polarity * flow_gain
    vol_raw = np.cumsum(flow_arr) * dt

    # 2. Capnogram calibration (mmHg)
    ch4_arr = None
    if ch4 is not None and len(ch4) == n:
        ch4_arr = (np.array(ch4) - co2_zero) * co2_gain

    # 3. Breath detection (end of expiration / start of inspiration)
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
                    if (time_arr[i] - time_arr[breath_starts[-1]]) > 0.8:
                        breath_starts.append(i)
                        in_expiration = False

    if not use_co2 or len(breath_starts) < 3:
        breath_starts = [0]
        in_neg = False
        for i in range(1, n):
            if flow_arr[i] < -20.0:
                in_neg = True
            elif in_neg and flow_arr[i] > 20.0:
                if (time_arr[i] - time_arr[breath_starts[-1]]) > 0.8:
                    breath_starts.append(i)
                    in_neg = False

    if breath_starts[-1] != n - 1:
        breath_starts.append(n - 1)

    # 4. Detrended Volume & Volumetric Capnography calculation
    vol_detrended = np.zeros_like(vol_raw)
    breaths = []

    for k in range(len(breath_starts) - 1):
        i0 = breath_starts[k]
        i1 = breath_starts[k + 1]
        cycle_flow = flow_arr[i0:i1 + 1]
        cycle_vol_raw = np.cumsum(cycle_flow) * dt
        drift = cycle_vol_raw[-1]
        correction = np.linspace(0, drift, len(cycle_vol_raw))
        cycle_vol_corr = cycle_vol_raw - correction
        
        if k == 0:
            vol_detrended[i0:i1 + 1] = cycle_vol_corr
        else:
            vol_detrended[i0 + 1:i1 + 1] = cycle_vol_corr[1:]

        peak_rel_idx = int(np.argmax(cycle_vol_corr))
        peak_idx = i0 + peak_rel_idx
        vt = float(cycle_vol_corr[peak_rel_idx])
        ti = float(time_arr[peak_idx] - time_arr[i0])
        te = float(time_arr[i1] - time_arr[peak_idx])
        tot_time = ti + te
        rr = float(60.0 / tot_time) if tot_time > 0 else 0.0

        et_co2 = 0.0
        fowler_vd = 0.0
        peco2 = 0.0
        paco2 = 0.0
        bohr_ratio = 0.0
        bohr_vd = 0.0
        vco2_ml = 0.0
        volcap_curve = []

        if ch4_arr is not None and vt > 50 and (i1 - peak_idx) > 5:
            co2_exp = ch4_arr[peak_idx:i1 + 1]
            v_exh = vt - cycle_vol_corr[peak_rel_idx:]
            et_co2 = float(np.max(co2_exp))
            
            sort_idx = np.argsort(v_exh)
            v_sorted = v_exh[sort_idx]
            co2_sorted = co2_exp[sort_idx]
            _, uniq_idx = np.unique(v_sorted, return_index=True)
            v_uniq = v_sorted[uniq_idx]
            co2_uniq = co2_sorted[uniq_idx]

            v_grid = np.linspace(0, vt, int(vt) + 1)
            co2_grid = np.interp(v_grid, v_uniq, co2_uniq)

            step_ds = max(1, len(v_grid) // 80)
            volcap_curve = [[round(float(v_grid[idx]), 1), round(float(co2_grid[idx]), 2)] 
                            for idx in range(0, len(v_grid), step_ds)]

            peco2 = float(np.mean(co2_grid))

            dco2_dv = np.gradient(co2_grid, v_grid)
            kernel_size = 11
            dco2_smooth = np.convolve(dco2_dv, np.ones(kernel_size) / kernel_size, mode='same')
            i_min = int(0.15 * len(v_grid))
            i_max = int(0.65 * len(v_grid))
            if i_max > i_min + 5:
                idx_inflect = i_min + int(np.argmax(dco2_smooth[i_min:i_max]))
                fowler_vd = float(v_grid[idx_inflect])
            else:
                fowler_vd = float(0.3 * vt)

            i_p3_start = min(int(fowler_vd + 0.15 * vt), int(0.65 * len(v_grid)))
            paco2 = float(np.mean(co2_grid[i_p3_start:])) if i_p3_start < len(co2_grid) else float(np.max(co2_grid))

            bohr_ratio = float((paco2 - peco2) / paco2) if paco2 > 0 else 0.0
            bohr_vd = float(bohr_ratio * vt)

            vco2_ml = float((peco2 / 760.0) * vt)

        breaths.append({
            "index": k + 1,
            "start_s": round(float(time_arr[i0]), 3),
            "peak_s": round(float(time_arr[peak_idx]), 3),
            "end_s": round(float(time_arr[i1]), 3),
            "start_idx": int(i0),
            "peak_idx": int(peak_idx),
            "end_idx": int(i1),
            "vt": round(vt, 1),
            "ti": round(ti, 3),
            "te": round(te, 3),
            "ie_ratio": f"1 : {round(te / ti, 2)}" if ti > 0 else "N/A",
            "rr": round(rr, 1),
            "et_co2": round(et_co2, 1),
            "fowler_vd": round(fowler_vd, 1),
            "peco2": round(peco2, 1),
            "paco2": round(paco2, 1),
            "bohr_ratio": round(bohr_ratio, 3),
            "bohr_vd": round(bohr_vd, 1),
            "vco2_ml": round(vco2_ml, 1),
            "volcap_curve": volcap_curve
        })

    return {
        "flow": flow_arr.tolist(),
        "ch4_mmhg": ch4_arr.tolist() if ch4_arr is not None else [],
        "vol_raw": vol_raw.tolist(),
        "vol_detrended": vol_detrended.tolist(),
        "breaths": breaths,
        "calibration": {
            "flow_gain": flow_gain,
            "co2_gain": co2_gain,
            "co2_zero": co2_zero,
            "baseline_offset": float(base)
        }
    }


class RespViewerHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(Path(__file__).resolve().parent), **kwargs)

    def do_GET(self):
        parsed = urlparse(self.path)
        path = parsed.path

        if path == "/api/files":
            self.send_json_response(self.get_available_files())
        elif path == "/api/load":
            query = parse_qs(parsed.query)
            filename = query.get("file", [None])[0]
            if not filename:
                self.send_error(400, "Missing 'file' parameter")
                return
            try:
                data = self.load_file_data(filename)
                self.send_json_response(data)
            except Exception as e:
                self.send_error(500, f"Error loading file: {str(e)}")
        elif path == "/api/process":
            query = parse_qs(parsed.query)
            filename = query.get("file", [None])[0]
            mode = query.get("mode", ["linear_detrend"])[0]
            polarity = float(query.get("polarity", [1.0])[0])
            flow_offset = query.get("offset", [None])[0]
            flow_gain = float(query.get("flow_gain", [19.56])[0])
            co2_zero = float(query.get("co2_zero", [22.0])[0])
            co2_gain = float(query.get("co2_gain", [0.231])[0])
            if flow_offset is not None:
                flow_offset = float(flow_offset)

            if not filename:
                self.send_error(400, "Missing 'file' parameter")
                return
            try:
                raw_data = self.load_file_data(filename)
                ch1 = raw_data["channels"].get("raw_ch1", [])
                ch4 = raw_data["channels"].get("raw_ch4", None)
                proc = compute_volume_and_breaths(
                    raw_data["time_s"],
                    ch1,
                    ch4=ch4,
                    mode=mode,
                    polarity=polarity,
                    flow_offset=flow_offset,
                    flow_gain=flow_gain,
                    co2_zero=co2_zero,
                    co2_gain=co2_gain,
                )
                self.send_json_response(proc)
            except Exception as e:
                self.send_error(500, f"Error processing file: {str(e)}")
        else:
            super().do_GET()

    def get_available_files(self):
        files = []
        ignored_dirs = {".git", "__pycache__", ".venv", "venv", "node_modules", "resp_viewer"}
        for p in sorted(WORKSPACE_DIR.rglob("*")):
            if any(part in ignored_dirs for part in p.parts):
                continue
            if p.is_file() and p.suffix.lower() in [".inp", ".csv"]:
                rel_path = p.relative_to(WORKSPACE_DIR).as_posix()
                size_kb = round(p.stat().st_size / 1024.0, 1)
                files.append({
                    "name": rel_path,
                    "type": p.suffix.lower()[1:],
                    "size_kb": size_kb,
                    "path": str(p.resolve()),
                })
        return {"files": files, "workspace": str(WORKSPACE_DIR)}

    def load_file_data(self, filename: str):
        target = (WORKSPACE_DIR / filename).resolve()
        if not str(target).startswith(str(WORKSPACE_DIR)) or not target.exists():
            raise FileNotFoundError(f"File not found: {filename}")
        if target.suffix.lower() == ".inp":
            data = decode_inp_file(target)
        elif target.suffix.lower() == ".csv":
            data = decode_csv_file(target)
        else:
            raise ValueError("Unsupported file extension.")

        ch1 = data["channels"].get("raw_ch1", [])
        ch4 = data["channels"].get("raw_ch4", None)
        proc = compute_volume_and_breaths(data["time_s"], ch1, ch4=ch4)
        data["computed"] = proc
        return data

    def send_json_response(self, obj):
        payload = json.dumps(obj).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(payload)

    def log_message(self, format, *args):
        sys.stderr.write(f"[{self.log_date_time_string()}] {format % args}\n")


def run_server(port=8088):
    server_address = ("0.0.0.0", port)
    httpd = HTTPServer(server_address, RespViewerHandler)
    print("=============================================================")
    print(f"  CapnoView szerver fut!")
    print(f"  Mappa: {WORKSPACE_DIR}")
    print(f"  Böngészőben megnyitható: http://localhost:{port}")
    print("=============================================================")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nSzerver leállítása.")
        httpd.server_close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CapnoView Server")
    parser.add_argument("-p", "--port", type=int, default=8088, help="Port to listen on")
    args = parser.parse_args()
    run_server(port=args.port)
