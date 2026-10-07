"""Verify a separately installed public package; write portable command reports."""
import argparse
import csv
import hashlib
import io
import json
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-tests", action="store_true",
                        help="Use when the same installed test suite has just been run separately")
    args = parser.parse_args()
    import nrcert
    origin = Path(nrcert.__file__).resolve()
    require(not origin.is_relative_to(ROOT / "src"),
            "Install the wheel before running this installed-package verification")
    logs = ROOT / "audit/commands"
    logs.mkdir(parents=True, exist_ok=True)
    output = ROOT / "output/verified-certificates"
    output.mkdir(parents=True, exist_ok=True)
    commands = []

    def run(name, arguments, expected=0):
        result = subprocess.run(arguments, cwd=ROOT, capture_output=True, text=True,
                                encoding="utf-8", errors="strict")
        (logs / f"{name}.txt").write_text(result.stdout + result.stderr, encoding="utf-8")
        commands.append({"name": name, "exit_code": result.returncode})
        require(result.returncode == expected,
                f"{name}: expected exit {expected}, got {result.returncode}: {result.stderr}")
        return result.stdout, result.stderr

    prefix = [sys.executable, "-I", "-B"]
    test_count = None
    if not args.skip_tests:
        _, test_log = run("installed-tests", prefix + ["-m", "unittest", "discover", "-s", "tests", "-v"])
        match = re.search(r"Ran (\d+) tests", test_log)
        require(match is not None and "\nOK\n" in test_log, "Test summary missing")
        test_count = int(match.group(1))
    csv_text, _ = run("comparison-reproduction", prefix + ["experiments/reproduce.py"])
    rows = list(csv.DictReader(io.StringIO(csv_text)))
    require(len(rows) == 30 and all(row["recomputed"] == "True" for row in rows),
            "Incomplete comparison reproduction")
    cli = Path(sys.executable).parent / ("nrcert.exe" if sys.platform == "win32" else "nrcert")
    require(cli.is_file(), "Installed nrcert executable is missing beside this Python interpreter")
    run("console-version", [str(cli), "--version"])
    for name in ("normal", "bidiagonal32", "coupled32", "shifted_j2",
                 "slow_contraction", "polynomial_segment"):
        status = 2 if name == "shifted_j2" else 0
        path = output / f"{name}.json"
        run(f"certify-{name}", [str(cli), "certify", str(ROOT / "examples" / f"{name}.json"),
                             "-o", str(path)], status)
        text, _ = run(f"check-{name}", [str(cli), "check", str(path)], status)
        require(json.loads(text)["valid"] is True, f"Check failed: {name}")
    report = {"status": "passed", "python": platform.python_version(),
              "platform": platform.platform(), "installed_import_verified": True,
              "tests_executed": not args.skip_tests, "tests_passed": test_count,
              "comparison_rows": len(rows), "console_examples_passed": 6,
              "package_source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                        for p in sorted(origin.parent.glob("*.py"))},
              "commands": commands}
    (ROOT / "audit/installed-package-verification.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
