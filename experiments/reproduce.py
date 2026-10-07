"""Print exact certificate comparisons as CSV; no output files are created.

Run after installing the package. Decimal q columns are display conversions,
not additional validated evidence. Computational timing is deliberately absent.
"""
import argparse
import csv
import json
import sys
from fractions import Fraction
from pathlib import Path
from nrcert import make_certificate, verify_certificate
from nrcert.reference import direct_power_reference

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--examples", type=Path,
                        default=Path(__file__).resolve().parents[1] / "examples")
    parser.add_argument("--direct-power-limit", type=int, default=64)
    args = parser.parse_args(argv)
    fields = [
        "fixture", "dimension", "resolution", "mode", "tolerance",
        "q_upper_squared", "q_decimal_display", "iterations",
        "one_step_induced_norm_iterations", "direct_power_iterations", "direct_power_status",
        "direct_power_limit", "status", "vertices", "recomputed",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for name in ("normal", "bidiagonal32", "coupled32", "shifted_j2", "slow_contraction"):
        request = json.loads((args.examples / f"{name}.json").read_text(encoding="utf-8"))
        direct = direct_power_reference(request, max_power=args.direct_power_limit)
        for resolution in (1, 4, 8):
            for mode in ("crouzeix_palencia", "established"):
                certificate = make_certificate({
                    **request, "direction_resolution": resolution, "mode": mode,
                })
                # Both Richardson comparisons are established theorems.
                checked = verify_certificate(certificate)
                evidence = certificate["evidence"]
                writer.writerow({
                    "fixture": name,
                    "dimension": len(request["matrix"]),
                    "resolution": resolution,
                    "mode": mode,
                    "tolerance": request["tolerance"],
                    "q_upper_squared": evidence["q_upper_squared"],
                    "q_decimal_display": float(Fraction(evidence["q_upper_squared"])) ** 0.5,
                    "iterations": evidence["iterations"],
                    "one_step_induced_norm_iterations": evidence["induced_norm_baseline"]["iterations"],
                    "direct_power_iterations": direct["first_certified_power"],
                    "direct_power_status": direct["status"],
                    "direct_power_limit": args.direct_power_limit,
                    "status": evidence["status"],
                    "vertices": len(certificate["vertices"]),
                    "recomputed": checked["valid"],
                })
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
