"""JSON command-line interface for exact numerical-range certificates."""
import argparse
import json
import sys
from pathlib import Path
from .core import CertificateError, make_certificate, verify_certificate

def _read(path):
    try:
        text = sys.stdin.read() if path == "-" else Path(path).read_text(encoding="utf-8")
        return json.loads(text)
    except (ValueError, UnicodeError, RecursionError) as exc:
        raise CertificateError("invalid or unsupported JSON input") from exc

def _write(value, path):
    text = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if path == "-":
        sys.stdout.write(text)
    else:
        Path(path).write_text(text, encoding="utf-8")

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="nrcert",
        description="Exact-rational polynomial and Richardson norm certificates.",
    )
    parser.add_argument("--version", action="version", version="nrcert 0.1.0")
    commands = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("certify", "construct a certificate from an exact-rational request"),
        ("check", "recompute every field of a stored certificate"),
    ):
        command = commands.add_parser(name, help=help_text)
        command.add_argument("input", help="JSON file, or - for standard input")
        command.add_argument("-o", "--output", default="-", help="JSON output file (default: stdout)")
        command.add_argument(
            "--allow-conditional", action="store_true",
            help="explicitly allow bounds conditional on the unestablished H2 assumption",
        )
    args = parser.parse_args(argv)
    try:
        value = _read(args.input)
        if args.command == "certify":
            if (isinstance(value, dict) and value.get("kind") == "polynomial"
                    and value.get("mode") == "conditional_two" and not args.allow_conditional):
                raise CertificateError("conditional H2 request requires --allow-conditional")
            certificate = make_certificate(value)
            result = verify_certificate(certificate, allow_conditional=args.allow_conditional)
            _write(certificate, args.output)
        else:
            result = verify_certificate(value, allow_conditional=args.allow_conditional)
            _write(result, args.output)
        return 0 if result["bound_certified"] else 2
    except (OSError, json.JSONDecodeError, CertificateError) as exc:
        print(f"nrcert: error: {exc}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())

