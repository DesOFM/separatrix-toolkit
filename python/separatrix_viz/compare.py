from __future__ import annotations

import numpy as np
from .io import load_table


def _xy(data: np.ndarray, schema: str) -> np.ndarray:
    if schema == "toolkit":
        if data.shape[1] < 6:
            raise ValueError("Toolkit schema requires at least 6 columns")
        return data[:, 4:6]
    if schema == "xy":
        if data.shape[1] < 2:
            raise ValueError("xy schema requires at least 2 columns")
        return data[:, 0:2]
    raise ValueError(f"Unknown schema: {schema}")


def directed_distances(a: np.ndarray, b: np.ndarray, chunk: int = 2048) -> np.ndarray:
    result = np.empty(len(a), dtype=float)
    for i in range(0, len(a), chunk):
        block = a[i:i+chunk]
        d2 = ((block[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)
        result[i:i+chunk] = np.sqrt(d2.min(axis=1))
    return result


def compare_files(reference: str, candidate: str, *, reference_schema: str = "xy",
                  candidate_schema: str = "toolkit") -> dict:
    ref = _xy(load_table(reference), reference_schema)
    cand = _xy(load_table(candidate), candidate_schema)
    d_rc = directed_distances(ref, cand)
    d_cr = directed_distances(cand, ref)
    return {
        "reference_to_candidate_rms": float(np.sqrt(np.mean(d_rc**2))),
        "candidate_to_reference_rms": float(np.sqrt(np.mean(d_cr**2))),
        "symmetric_hausdorff": float(max(d_rc.max(), d_cr.max())),
    }
