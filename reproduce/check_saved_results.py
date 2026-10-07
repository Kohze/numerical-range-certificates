"""Recompute retained exact data without changing files or remeasuring timings."""
import csv
import importlib.util
import io
import json
import statistics
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    from nrcert import verify_certificate
    run = subprocess.run([sys.executable, "-I", "-B", str(ROOT / "experiments/reproduce.py")],
                         cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
    require(run.returncode == 0, run.stderr)
    actual = list(csv.DictReader(io.StringIO(run.stdout)))
    expected = csv_rows(ROOT / "results/comparisons.csv")
    require(len(actual) == 30 and actual == expected, "Comparison CSV differs from recomputation")
    checked = 0
    for path in sorted((ROOT / "results/certificates").glob("*.json")):
        certificate = json.loads(path.read_text(encoding="utf-8"))
        require(verify_certificate(certificate)["valid"], f"Invalid saved certificate: {path.name}")
        checked += 1
    require(checked == 6, "Expected six saved certificates")
    module = load_module("public_ellipse_study", "experiments/ellipse_error_budget.py")
    calculated = module.study()
    saved = csv_rows(ROOT / "results/ellipse-error-budget.csv")
    require(len(calculated) == len(saved) == 144, "Ellipse row count differs")
    for number, (row, old) in enumerate(zip(calculated, saved), 1):
        require(set(row) == set(old), f"Ellipse schema differs at row {number}")
        for key, value in row.items():
            if key.endswith("_display"):
                continue
            require(str(value) == old[key], f"Exact ellipse field differs at row {number}: {key}")
    costs = json.loads((ROOT / "results/certificate-cost-study.json").read_text(encoding="utf-8"))
    cost_rows = csv_rows(ROOT / "results/certificate-cost-study.csv")
    require(len(cost_rows) == len(costs["observations"]) == 10, "Timing fixture count differs")
    require(costs["repetitions"] == 3, "Unexpected recorded timing repetitions")
    medians = 0
    for row, observation in zip(cost_rows, costs["observations"]):
        require(row["fixture"] == observation["fixture"], "Timing fixture order differs")
        for kind in ("construction", "checking", "norm_only", "direct_power"):
            times = observation[kind + "_seconds"]
            require(len(times) == 3 and all(value >= 0 for value in times),
                    "Invalid recorded timing observations")
            require(float(row[kind + "_median_seconds"]) == statistics.median(times),
                    f"Recorded median differs: {row['fixture']}/{kind}")
            medians += 1
    print(json.dumps({"status": "passed", "comparison_rows": len(actual),
                      "certificates_checked": checked, "exact_ellipse_rows": len(saved),
                      "recorded_timing_medians_checked": medians,
                      "timings_remeasured": False}))


if __name__ == "__main__":
    main()
