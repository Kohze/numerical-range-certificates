"""Print an exact direct-power reference for a Richardson JSON request."""
import argparse
import json
from pathlib import Path
from nrcert.reference import direct_power_reference

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--max-power", type=int, default=64)
    args = parser.parse_args(argv)
    request = json.loads(args.input.read_text(encoding="utf-8"))
    print(json.dumps(direct_power_reference(request, args.max_power), indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
