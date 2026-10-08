#!/usr/bin/env python3
"""Grid-convergence study for Separatrix Toolkit.

Generates two analytic benchmarks with known separatrices, runs a compiled
Separatrix Toolkit executable, and evaluates saddle, level-set, and geometric
convergence metrics.

Requires: numpy, shapely.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import time
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, MultiLineString
from shapely import points, distance


def duffing_psi(x, y, x0, y0):
    X = x - x0
    Y = y - y0
    return 0.5 * Y**2 - 0.5 * X**2 + 0.25 * X**4


def heteroclinic_psi(x, y, x0, y0):
    X = x - x0
    Y = y - y0
    return (X**2 - 1.0)**2 - Y**2


CASES = {
    "duffing": {
        "domain": (-2.0, 2.0, -1.2, 1.2),
        "shift": (0.137, -0.083),
        "psi": duffing_psi,
        "exact_saddles": lambda x0, y0: np.array([[x0, y0]]),
    },
    "heteroclinic": {
        "domain": (-2.0, 2.0, -4.0, 4.0),
        "shift": (0.137, -0.083),
        "psi": heteroclinic_psi,
        "exact_saddles": lambda x0, y0: np.array([[x0 - 1.0, y0], [x0 + 1.0, y0]]),
    },
}


def exact_geometry_and_points(name: str, n: int = 20000):
    case = CASES[name]
    x0, y0 = case["shift"]

    if name == "duffing":
        X = np.linspace(-math.sqrt(2.0), math.sqrt(2.0), n)
        Y = np.sqrt(np.maximum(0.0, X**2 - 0.5 * X**4))
        p1 = np.c_[X + x0, Y + y0]
        p2 = np.c_[X + x0, -Y + y0]
    else:
        xmin, xmax, _, _ = case["domain"]
        x = np.linspace(xmin, xmax, n)
        X = x - x0
        Y = X**2 - 1.0
        p1 = np.c_[x, Y + y0]
        p2 = np.c_[x, -Y + y0]

    pts = np.vstack([p1, p2])
    geom = MultiLineString([LineString(p1), LineString(p2)])
    return geom, pts


def read_table(path):
    a = np.loadtxt(path, comments="#")
    return a.reshape(1, -1) if a.ndim == 1 else a


def candidate_geometry(table):
    lines = []
    for sep_id in np.unique(table[:, 1].astype(int)):
        block = table[table[:, 1].astype(int) == sep_id]
        for branch_id in np.unique(block[:, 2].astype(int)):
            branch = block[block[:, 2].astype(int) == branch_id]
            branch = branch[np.argsort(branch[:, 3])]
            if len(branch) >= 2:
                lines.append(LineString(branch[:, 4:6]))
    return MultiLineString(lines)


def write_field(path, name, n):
    case = CASES[name]
    xmin, xmax, ymin, ymax = case["domain"]
    x0, y0 = case["shift"]
    x = np.linspace(xmin, xmax, n)
    y = np.linspace(ymin, ymax, n)

    with Path(path).open("w") as f:
        for xv in x:
            values = case["psi"](xv, y, x0, y0)
            for yv, value in zip(y, values):
                f.write(f"{xv:.16e} {yv:.16e} {value:.16e}\n")
    return x, y


def saddle_error(found, exact):
    remaining = set(range(len(found)))
    errors = []
    for target in exact:
        d = np.linalg.norm(found[:, 2:4] - target, axis=1)
        for idx in np.argsort(d):
            if int(idx) in remaining:
                remaining.remove(int(idx))
                errors.append(float(d[idx]))
                break
    return max(errors)


def farfield_metrics(candidate_geom, candidate_points, exact_geom, exact_points,
                     saddles, radius=0.2):
    """Sampled symmetric geometric error outside fixed saddle disks.

    Candidate vertices are measured against the dense analytic geometry, while
    dense analytic reference points are measured against the candidate polyline.
    With a sufficiently dense analytic sampling this is an accurate Hausdorff
    estimate for the convergence study and is much faster than repeated exact
    polygonal Hausdorff densification.
    """
    cdist = np.min(
        np.linalg.norm(candidate_points[:, None, :] - saddles[None, :, :], axis=2),
        axis=1,
    )
    edist = np.min(
        np.linalg.norm(exact_points[:, None, :] - saddles[None, :, :], axis=2),
        axis=1,
    )
    cp = candidate_points[cdist > radius]
    ep = exact_points[edist > radius]

    dc = np.asarray(distance(points(cp[:, 0], cp[:, 1]), exact_geom))
    de = np.asarray(distance(points(ep[:, 0], ep[:, 1]), candidate_geom))

    hd = max(float(dc.max()), float(de.max()))
    rms = math.sqrt(0.5 * (float(np.mean(dc**2)) + float(np.mean(de**2))))
    return hd, rms


def fitted_order(h, e):
    return float(np.polyfit(np.log(h), np.log(e), 1)[0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", required=True, help="compiled separatrix_toolkit executable")
    parser.add_argument("--output", default="convergence_runs")
    parser.add_argument("--sizes", nargs="+", type=int, default=[101, 201, 401, 801])
    parser.add_argument("--far-radius", type=float, default=0.2)
    parser.add_argument("--reference-points", type=int, default=80001,
                        help="points per analytic branch")
    args = parser.parse_args()

    outroot = Path(args.output)
    outroot.mkdir(parents=True, exist_ok=True)
    all_rows = []

    for name, case in CASES.items():
        exact_geom, exact_pts = exact_geometry_and_points(name, args.reference_points)
        exact_saddles = case["exact_saddles"](*case["shift"])

        for n in args.sizes:
            work = outroot / name / f"N{n}"
            work.mkdir(parents=True, exist_ok=True)
            x, y = write_field(work / "field.dat", name, n)

            cmd = [
                args.solver, str(work / "field.dat"),
                "--saddles-output", str(work / "saddles.dat"),
                "--separatrices-output", str(work / "separatrices.dat"),
                "--diagnostics-output", str(work / "diagnostics.txt"),
                "--quiet",
            ]
            t0 = time.perf_counter()
            subprocess.run(cmd, check=True)
            runtime = time.perf_counter() - t0

            saddles = read_table(work / "saddles.dat")
            sep = read_table(work / "separatrices.dat")
            cand_pts = sep[:, 4:6]
            cand_geom = candidate_geometry(sep)

            x0, y0 = case["shift"]
            level_error = np.abs(case["psi"](cand_pts[:, 0], cand_pts[:, 1], x0, y0))
            far_hd, far_rms = farfield_metrics(
                cand_geom, cand_pts, exact_geom, exact_pts, exact_saddles,
                radius=args.far_radius,
            )

            all_rows.append({
                "case": name,
                "N": n,
                "hx": float(x[1] - x[0]),
                "hy": float(y[1] - y[0]),
                "h": float(math.hypot(x[1] - x[0], y[1] - y[0])),
                "saddle_max_error": saddle_error(saddles, exact_saddles),
                "exact_level_rms": float(np.sqrt(np.mean(level_error**2))),
                "exact_level_max": float(level_error.max()),
                "farfield_hausdorff_estimate": far_hd,
                "farfield_rms": far_rms,
                "runtime_s": runtime,
            })

    for name in CASES:
        block = [r for r in all_rows if r["case"] == name]
        h = [r["h"] for r in block]
        for metric in [
            "saddle_max_error", "exact_level_rms", "exact_level_max",
            "farfield_hausdorff_estimate", "farfield_rms",
        ]:
            p = fitted_order(h, [r[metric] for r in block])
            for r in block:
                r[f"fit_order_{metric}"] = p

    csv_path = outroot / "convergence_results.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)

    print(json.dumps(all_rows, indent=2))
    print(f"Saved: {csv_path}")


if __name__ == "__main__":
    main()
