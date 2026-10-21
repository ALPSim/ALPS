# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
"""Plot fresh nitrogen/scaling JSON without running a calculation.

Replaces the backend-specific report generators in skilledwolf/alpsx at
1adba871e plus local edits. Uses their saved-data workflow, with no historical
numbers, fixed bond schedules or inferred convergence claims.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, help="defaults to the input name with .png")
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    rows = data["dmrg"]
    if not rows:
        parser.error("the input has no completed DMRG runs")
    bonds = [r["bond_requested"] for r in rows]
    nitrogen = "fidelity_vs_exact" in rows[0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    if nitrogen:
        axes[0].plot(bonds, [r["error_vs_exact_mEh"] for r in rows], "o-", color="#087e8b")
        axes[0].set(ylabel="Energy error / mEh", title="Final-state energy versus FCI")
        axes[1].plot(bonds, [r["fidelity_vs_exact"] for r in rows], "o-", color="#c76d24")
        axes[1].set(ylabel="Squared overlap with FCI", ylim=(0, 1.02), title="State fidelity")
    else:
        axes[0].plot(bonds, [r["energy"] for r in rows], "o-", color="#087e8b")
        axes[0].set(ylabel=f"Energy / {data['problem']['energy_unit']}", title=data["problem"]["name"])
        axes[1].plot(bonds, [r["seconds_including_rdms"] for r in rows], "o-", color="#c76d24")
        axes[1].set(ylabel="Seconds including RDMs", title="Single-run solver time")
    for ax in axes:
        ax.set_xlabel("Requested retained states")
        ax.spines[["top", "right"]].set_visible(False)
    output = args.output or args.input.with_suffix(".png")
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
