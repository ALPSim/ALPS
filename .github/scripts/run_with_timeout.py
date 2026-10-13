"""Stop native test processes even when extension code holds Python's GIL."""

import argparse
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timeout", type=float, help="maximum runtime in seconds")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.command or args.timeout <= 0:
        parser.error("a positive timeout and a command are required")
    try:
        return subprocess.run(args.command, timeout=args.timeout).returncode
    except subprocess.TimeoutExpired:
        print(f"Command exceeded its {args.timeout:g}-second timeout", file=sys.stderr)
        return 124


if __name__ == "__main__":
    raise SystemExit(main())
