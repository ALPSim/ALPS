"""Route the simulation-help form using its exact Method choices."""

import json
import os
from pathlib import Path


def assignee(body):
    routes = {
        "Exact Diagonalization (ED)": "vws100",
        "Quantum Monte Carlo (SSE)": "wistaria",
        "Quantum Monte Carlo (Worm)": "LodePollet",
        "Density Matrix Renormalization Group (DMRG)": "afeiguin",
        "Dynamical Mean Field Theory (DMFT)": "egull",
    }
    in_method = False
    for line in (body or "").splitlines():
        line = line.strip()
        if line.startswith("### "):
            in_method = line == "### Method"
        elif in_method and line:
            return routes.get(line, "Ooolab")
    return "Ooolab"


if __name__ == "__main__":
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text(encoding="utf-8"))
    with Path(os.environ["GITHUB_OUTPUT"]).open("a", encoding="utf-8") as output:
        output.write(f"assignee={assignee(event['issue'].get('body'))}\n")
