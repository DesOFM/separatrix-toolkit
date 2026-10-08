from __future__ import annotations

from pathlib import Path
import numpy as np


def load_table(path: str | Path) -> np.ndarray:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    data = np.loadtxt(p, comments="#")
    if data.size == 0:
        return np.empty((0, 0), dtype=float)
    return data.reshape(1, -1) if data.ndim == 1 else data


def load_optional(path: str | Path | None) -> np.ndarray | None:
    if path is None:
        return None
    p = Path(path)
    if not p.exists():
        return None
    data = load_table(p)
    return None if data.size == 0 else data


def prepare_field(source: str | Path, output: str | Path, *, x_col: int = 1,
                  y_col: int = 2, psi_col: int = 3) -> None:
    data = load_table(source)
    cols = [x_col - 1, y_col - 1, psi_col - 1]
    if min(cols) < 0 or max(cols) >= data.shape[1]:
        raise ValueError(f"Requested columns exceed source width ({data.shape[1]})")
    out = data[:, cols]
    np.savetxt(output, out, fmt="%.15e", header="x y psi")


def separatrix_summary(sep: np.ndarray) -> dict:
    if sep.size == 0:
        return {"layers": 0, "networks": 0, "branches": 0, "points": 0}
    layers = np.unique(sep[:, 0].astype(int))
    networks = {(int(row[0]), int(row[1])) for row in sep}
    branches = {(int(row[0]), int(row[1]), int(row[2])) for row in sep}
    return {
        "layers": len(layers),
        "networks": len(networks),
        "branches": len(branches),
        "points": len(sep),
    }
