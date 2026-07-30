#!/usr/bin/env python3
"""Decode the reverse-engineered Winlung-style .inp recording format."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import numpy as np
import pandas as pd


def read_pascal_u16_string(blob: bytes, pos: int) -> tuple[str, int]:
    if pos + 2 > len(blob):
        raise ValueError("Truncated string length.")
    length = struct.unpack_from("<H", blob, pos)[0]
    pos += 2
    if pos + length > len(blob):
        raise ValueError("Truncated string content.")
    value = blob[pos:pos + length].decode("ascii", errors="replace")
    return value, pos + length


def decode_inp(path: Path) -> tuple[dict, np.ndarray]:
    blob = path.read_bytes()
    if len(blob) < 44:
        raise ValueError("File is too short.")

    header = {
        "sample_rate_hz": struct.unpack_from("<f", blob, 0)[0],
        "sample_count": struct.unpack_from("<H", blob, 4)[0],
        "channel_count": struct.unpack_from("<H", blob, 6)[0],
        "duration_s": struct.unpack_from("<f", blob, 8)[0],
        "reserved_24": blob[12:36],
        "marker": struct.unpack_from("<i", blob, 36)[0],
        "flags": struct.unpack_from("<H", blob, 40)[0],
    }

    pos = 42
    header["recorded_at"], pos = read_pascal_u16_string(blob, pos)
    header["subject_id"], pos = read_pascal_u16_string(blob, pos)

    if pos + 2 > len(blob):
        raise ValueError("Truncated metadata.")
    header["unknown_u16"] = struct.unpack_from("<H", blob, pos)[0]
    pos += 2

    header["signal_configuration_path"], pos = read_pascal_u16_string(blob, pos)
    header["data_offset"] = pos

    payload = np.frombuffer(blob, dtype="<i2", offset=pos)
    expected = header["sample_count"] * header["channel_count"]
    if payload.size != expected:
        raise ValueError(
            f"Payload has {payload.size} int16 values; expected {expected}."
        )

    signals = payload.reshape(
        header["sample_count"], header["channel_count"]
    )
    return header, signals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=None)
    args = parser.parse_args()

    header, signals = decode_inp(args.input_file)
    output = args.output or args.input_file.with_suffix(".csv")

    time_s = np.arange(header["sample_count"]) / header["sample_rate_hz"]
    df = pd.DataFrame({"time_s": time_s})
    for index in range(header["channel_count"]):
        df[f"raw_ch{index + 1}"] = signals[:, index]
    df.to_csv(output, index=False)

    print("Header:")
    for key, value in header.items():
        if key != "reserved_24":
            print(f"  {key}: {value}")
    print(f"Wrote: {output}")


if __name__ == "__main__":
    main()
