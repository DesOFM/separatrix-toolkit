from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from .io import load_table, load_optional


def _field_grid(field: np.ndarray):
    x = np.unique(field[:, 0])
    y = np.unique(field[:, 1])
    if len(x) * len(y) != len(field):
        return None
    z = field[:, 2].reshape(len(x), len(y))
    return x, y, z.T


def plot_networks(
    separatrices: str | Path,
    *,
    saddles: str | Path | None = None,
    field: str | Path | None = None,
    layer: int | None = None,
    title: str | None = None,
    output: str | Path | None = None,
    formats: list[str] | None = None,
    dpi: int = 300,
    contour_levels: int = 25,
    annotate_saddles: bool = True,
    equal_aspect: bool = True,
    show: bool = False,
) -> list[Path]:
    sep = load_table(separatrices)
    if sep.size == 0:
        raise ValueError("No separatrix data found")

    saddle_data = load_optional(saddles)
    field_data = load_optional(field)

    if layer is not None:
        sep = sep[sep[:, 0].astype(int) == layer]
        if saddle_data is not None:
            saddle_data = saddle_data[saddle_data[:, 0].astype(int) == layer]
        if sep.size == 0:
            raise ValueError(f"No separatrix data for layer {layer}")

    fig, ax = plt.subplots(figsize=(10, 7))

    if field_data is not None and field_data.shape[1] >= 3:
        grid = _field_grid(field_data)
        if grid is not None:
            xg, yg, zg = grid
            ax.contour(xg, yg, zg, levels=contour_levels, linewidths=0.55, alpha=0.55)

    layer_ids = np.unique(sep[:, 0].astype(int))
    for current_layer in layer_ids:
        lblock = sep[sep[:, 0].astype(int) == current_layer]
        for sid in np.unique(lblock[:, 1].astype(int)):
            block = lblock[lblock[:, 1].astype(int) == sid]
            first = True
            for bid in np.unique(block[:, 2].astype(int)):
                branch = block[block[:, 2].astype(int) == bid]
                branch = branch[np.argsort(branch[:, 3])]
                label = f"L{current_layer} / separatrix {sid}" if first else None
                ax.plot(branch[:, 4], branch[:, 5], linewidth=2.0, label=label)
                first = False

    if saddle_data is not None and saddle_data.size > 0 and saddle_data.shape[1] >= 4:
        ax.scatter(saddle_data[:, 2], saddle_data[:, 3], zorder=5, label="Saddles")
        if annotate_saddles:
            for row in saddle_data:
                ax.annotate(
                    f"S{int(row[1])}",
                    (row[2], row[3]),
                    xytext=(5, 5),
                    textcoords="offset points",
                )

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title or "Separatrix network")
    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    outputs: list[Path] = []
    if output is not None:
        out = Path(output)
        if formats:
            stem = out.with_suffix("")
            for fmt in formats:
                target = stem.with_suffix("." + fmt.lstrip("."))
                fig.savefig(target, dpi=dpi, bbox_inches="tight")
                outputs.append(target)
        else:
            fig.savefig(out, dpi=dpi, bbox_inches="tight")
            outputs.append(out)

    if show:
        plt.show()
    else:
        plt.close(fig)

    return outputs
