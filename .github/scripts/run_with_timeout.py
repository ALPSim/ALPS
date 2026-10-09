"""Stop native test processes even when extension code holds Python's GIL."""

import argparse
import os
import signal
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timeout", type=float, help="maximum runtime in seconds")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.command or args.timeout <= 0:
        parser.error("a positive timeout and a command are required")
    # CI runs on POSIX. Keep nested runners in a group we can terminate together.
    with subprocess.Popen(args.command, start_new_session=True) as process:
        try:
            return process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            try:
                # MPI launchers need a chance to stop ranks in separate groups.
                os.killpg(process.pid, signal.SIGTERM)
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
            except ProcessLookupError:
                pass
            # Also stop children that outlive their launcher or ignore SIGTERM.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            print(f"Command exceeded its {args.timeout:g}-second timeout", file=sys.stderr)
            return 124


if __name__ == "__main__":
    raise SystemExit(main())
