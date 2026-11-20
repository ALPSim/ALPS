# SPDX-License-Identifier: MIT
"""Hypothesis tests of Monte Carlo averages against exact values.

This follows Wallerberger and Gull, Phys. Rev. E 96, 053303 (2017),
https://doi.org/10.1103/PhysRevE.96.053303, Sec. II A. Each average is
compared with its exact value through z = (mean - exact) / error, where the
error is the binning error that the application reports. Many bins enter that
error, so z is taken as standard normal and the two-sided p value is
erfc(|z| / sqrt(2)).

All comparisons of one test form a family that shares the level ALPHA: each of
the K comparisons is rejected when p < ALPHA / K (Bonferroni). The runs are
seeded, so a test reproduces one trajectory; ALPHA is the probability that a
correct code fails the test after a change that alters the random stream.

A bias of the size of the error bar is not caught by one run. The tests choose
the number of sweeps so that the bias they are meant to catch is several
error bars.
"""
import math
import xml.etree.ElementTree as ET

ALPHA = 1e-3


def read_averages(path):
    """Return {name: (mean, error, converged)} from the task's AVERAGES in a task output file.

    The averages of the individual runs (under MCRUN) are not read.
    """
    averages = {}
    task = ET.parse(path).getroot().find("AVERAGES")
    for average in task.findall("SCALAR_AVERAGE") if task is not None else []:
        error = average.find("ERROR")
        averages[average.get("name")] = (float(average.findtext("MEAN")), float(error.text),
                                         error.get("converged", "yes"))
    return averages


def compare(label, averages, exact, alpha=ALPHA, family=None):
    """Test each average named in `exact` against its exact value.

    Prints one line per comparison and returns a list of failure messages. An
    average is a failure if it is missing, if its error is not a positive
    finite number, if the binning analysis reports it as not converged, or if
    its two-sided p value is below alpha / family. `family` is the number of
    comparisons sharing alpha; it defaults to len(exact).
    """
    family = family or len(exact)
    threshold = alpha / family
    failures = []
    for name, value in exact.items():
        if name not in averages:
            failures.append(f"{label}: {name} missing from the output")
            continue
        mean, error, converged = averages[name]
        if not (math.isfinite(error) and error > 0):
            failures.append(f"{label}: {name} has error {error}")
            continue
        z = (mean - value) / error
        p = math.erfc(abs(z) / math.sqrt(2))
        note = "" if converged == "yes" else f", error convergence: {converged}"
        print(f"{label}: {name} = {mean:.6g} +- {error:.3g}, exact {value:.6g}, z = {z:+.2f}, p = {p:.2g}{note}")
        if converged == "no":
            failures.append(f"{label}: {name} error not converged")
        if p < threshold:
            failures.append(f"{label}: {name} = {mean:.6g} +- {error:.3g} differs from exact {value:.6g} "
                            f"by {z:+.2f} error bars (p = {p:.2g} < {threshold:.2g})")
    return failures
