"""Run required MPI tests with a distinct JUnit report and exit code per rank."""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("tests", nargs="+")
    args = parser.parse_args()
    # These are required imports: a promised MPI job must never silently skip.
    from mpi4py import MPI
    import pytest

    world = MPI.COMM_WORLD
    if world.size != args.ranks or args.ranks < 2:
        raise RuntimeError(f"Expected {args.ranks} MPI ranks (at least two), got {world.size}")
    args.output.mkdir(parents=True, exist_ok=True)
    result = pytest.main([
        "-q", *args.tests,
        f"--junitxml={args.output / f'rank-{world.rank}.xml'}",
    ])
    return world.allreduce(int(result), op=MPI.MAX)


if __name__ == "__main__":
    raise SystemExit(main())
