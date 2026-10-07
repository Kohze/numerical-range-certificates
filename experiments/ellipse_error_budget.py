"""Exact finite-error study on an analytically known nonnormal numerical range.

The experimental analytic support producer is not a shipped backend. Geometry
helpers are shared with nrcert; this is not an independent software verification.
All assertions compare exact Fractions. Decimal CSV columns are display only.
"""
import argparse
import csv
import json
from fractions import Fraction as F
from pathlib import Path

from nrcert import sqrt_upper
from nrcert.core import _directions, _supports, _vertices, _edges, _evaluate, _norm_squared


def study():
    matrix = [[F(-1), F(1)], [F(0), F(1)]]
    rows = []
    for bits in (8, 24, 40):
        step = F(1, 2**bits)
        radius = sqrt_upper(F(5, 4), bits)
        for resolution in (1, 2, 4, 8):
            for producer in ("analytic_ellipse", "gershgorin"):
                if producer == "gershgorin":
                    supports = _supports(matrix, resolution, bits)
                else:
                    supports = [(a, b, sqrt_upper(F(5*a*a + b*b, 4), bits), [])
                                for a, b in _directions(resolution)]
                eta = max((upper - max(F(0), sqrt_upper(F(5*a*a+b*b, 4), bits)-step)) / resolution
                          for a, b, upper, _ in supports)
                # ||(a,b)||>=resolution. Delta<=atan(1/k)<=1/k. Rational
                # cos/sin bounds replace numerical trigonometric calculations.
                rho = (radius + eta) / (1-F(1, 8*resolution**2))
                excess = eta + (radius+rho)/(2*resolution)
                vertices = _vertices(supports)
                if not all(_norm_squared(v) <= rho*rho for v in vertices):
                    raise ArithmeticError("a polygon vertex exceeds the proved radius")
                edge_squared = max(_norm_squared((b[0]-a[0], b[1]-a[1])) for a,b in _edges(vertices))
                for polynomial, coefficients in (("z_squared", [F(0),F(0),F(1)]),
                                                  ("one_minus_z_squared", [F(1),F(0),F(-1)])):
                    derivative = 2*rho
                    for segments in (4,16,64):
                        maximum_squared = F(0)
                        samples = 0
                        for first,last in _edges(vertices):
                            for j in range(segments+1):
                                z = tuple(first[k]+(last[k]-first[k])*F(j,segments) for k in (0,1))
                                maximum_squared = max(maximum_squared, _norm_squared(_evaluate(coefficients,z)))
                                samples += 1
                        sample_upper = sqrt_upper(maximum_squared,bits)
                        delta = sqrt_upper(edge_squared,bits)/(2*segments)
                        mesh_term = delta*derivative
                        upper = sample_upper+mesh_term
                        known_supremum = F(5,4)
                        observed_excess = upper-known_supremum
                        budget = derivative*(excess+delta)+step
                        if not 0 <= observed_excess <= budget:
                            raise ArithmeticError("the exact scalar error-budget inequality failed")
                        rational_values = {"support_slack_upper":eta,"radius_upper":rho,
                            "hausdorff_excess_upper":excess,"derivative_upper":derivative,
                            "mesh_radius_upper":delta,"sample_upper":sample_upper,
                            "mesh_term":mesh_term,"U_Q":upper,"M_Q":known_supremum,
                            "observed_excess":observed_excess,"excess_budget":budget}
                        row = {"producer":producer,"polynomial":polynomial,"sqrt_bits":bits,
                               "resolution":resolution,"segments_per_edge":segments,
                               "vertices":len(vertices),"sample_count":samples}
                        row.update({key:str(value) for key,value in rational_values.items()})
                        row.update({key+"_display":float(value) for key,value in rational_values.items()})
                        rows.append(row)
    return rows


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().parents[1]/"results/ellipse-error-budget.csv")
    args=parser.parse_args()
    rows=study()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary={"rows":len(rows),"all_exact_budget_inequalities_passed":True,
             "matrix":[[-1,1],[0,1]],"numerical_range":"x^2/(5/4)+y^2/(1/4)<=1",
             "known_scalar_supremum":"5/4 for both polynomials",
             "exact_matrix_norms":{"z_squared":"1","one_minus_z_squared":"0"},
             "analytic_support_producer_shipped":False,"independent_implementation":False,
             "study_scope":"finite scalar diagnostic; not a performance or universal theorem proof"}
    args.output.with_suffix(".json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
