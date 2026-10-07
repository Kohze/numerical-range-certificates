"""Small measured study of certificate work; it does not time linear solvers.

Reports all three observations, their median, and serialized evidence sizes.
The direct-power run has a fixed cap and may fail to certify. Such runs must
not be ranked as successful certificates. Elapsed timings are environment-specific.
"""
import argparse
import csv
import hashlib
import json
import platform
import statistics
import time
from fractions import Fraction as F
from pathlib import Path

from nrcert import make_certificate, verify_certificate
from nrcert.reference import direct_power_reference
from nrcert.core import _normalise, _induced_product, _minimum_iterations

ROOT=Path(__file__).resolve().parents[1]


def coupled(n):
    return [["1" if i==j and i%2==0 else "11/10" if i==j else
             "-4/5" if j==i+1 and i%2==0 else "-1/100" if j==i+1 else "0"
             for j in range(n)] for i in range(n)]


def discretization(n):
    # -nu*u'' + velocity*u' + reaction*u on (0,1), centered differences,
    # homogeneous Dirichlet boundaries, fixed nu=1/100, velocity=1/10,
    # reaction=1, and h=1/(n+1). This is a stated model, not a dataset.
    h=F(1,n+1)
    diffusion=F(1,100)/(h*h)
    advection=F(1,10)/(2*h)
    diagonal=1+2*diffusion
    A=[[str(diagonal) if i==j else str(-diffusion-advection) if j==i-1 else
        str(-diffusion+advection) if j==i+1 else "0" for j in range(n)] for i in range(n)]
    return A,str(1/diagonal)


def requests():
    for name in ("normal","bidiagonal32","coupled32","shifted_j2","slow_contraction"):
        request=json.loads((ROOT/"examples"/(name+".json")).read_text(encoding="utf-8"))
        yield name,request,"analytic control"
    yield "nilpotent2_099",{"matrix":[[1,"99/100"],[0,1]],"tolerance":"1/100000000"},"analytic counterexample control"
    yield "coupled8",{"matrix":coupled(8),"tolerance":"1/100000000"},"analytic control"
    for n in (8,16,32):
        A,alpha=discretization(n)
        yield f"reaction_advection_diffusion{n}",{"matrix":A,"alpha":alpha,"tolerance":"1/100000000"},"specified finite-difference model"


def measured(function,repeats):
    times=[]
    result=None
    for _ in range(repeats):
        start=time.perf_counter_ns()
        result=function()
        times.append((time.perf_counter_ns()-start)/1e9)
    return result,times


def norm_only_reference(raw):
    """Time the existing one-step criterion without constructing supports."""
    request,A=_normalise(raw)
    n=len(A)
    alpha=F(request["alpha"])
    B=[[F(int(i==j))-alpha*A[i][j] for j in range(n)] for i in range(n)]
    squared=_induced_product(B)
    iterations,status=_minimum_iterations(squared,F(1),F(request["tolerance"])**2,request["max_iterations"])
    return {"iterations":iterations,"status":status,"norm_product":str(squared)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=ROOT/"results/certificate-cost-study.csv")
    parser.add_argument("--repeats",type=int,default=3)
    args=parser.parse_args()
    if args.repeats<3:
        parser.error("at least three observations are required")
    make_certificate({"matrix":[[1]]})  # one untimed import/interpreter warmup
    rows=[]
    commands=[]
    for name,raw,category in requests():
        request={**raw,"mode":"established","direction_resolution":1,"sqrt_bits":40}
        certificate,construction_times=measured(lambda:make_certificate(request),args.repeats)
        checked,checking_times=measured(lambda:verify_certificate(certificate),args.repeats)
        norm_only,norm_times=measured(lambda:norm_only_reference(request),args.repeats)
        direct,power_times=measured(lambda:direct_power_reference(request,max_power=64),args.repeats)
        if not checked["valid"] or norm_only["iterations"] != certificate["evidence"]["induced_norm_baseline"]["iterations"]:
            raise ArithmeticError("certificate or one-step reference comparison failed")
        e=certificate["evidence"]
        row={"fixture":name,"category":category,"dimension":len(request["matrix"]),
             "mode":request["mode"],"numerical_range_iterations":e["iterations"],
             "one_step_iterations":e["induced_norm_baseline"]["iterations"],
             "direct_power_iterations":direct["first_certified_power"],"direct_power_status":direct["status"],
             "direct_power_cap":64,"construction_median_seconds":statistics.median(construction_times),
             "checking_median_seconds":statistics.median(checking_times),
             "norm_only_median_seconds":statistics.median(norm_times),
             "direct_power_median_seconds":statistics.median(power_times),
             "certificate_json_bytes":len(json.dumps(certificate,sort_keys=True,separators=(",",":"),allow_nan=False).encode()),
             "direct_power_json_bytes":len(json.dumps(direct,sort_keys=True,separators=(",",":"),allow_nan=False).encode())}
        rows.append(row)
        commands.append({"fixture":name,"request":request,"construction_seconds":construction_times,
                         "checking_seconds":checking_times,"direct_power_seconds":power_times,
                         "norm_only_seconds":norm_times,"norm_only_result":norm_only,
                         "assumption":certificate["assumption"],"outcomes":row})
        print(f"Completed {name}: NR={e['iterations']}, direct={direct['first_certified_power']} ({direct['status']})",flush=True)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report={"python":platform.python_version(),"platform":platform.platform(),"processor":platform.processor(),
            "repetitions":args.repeats,"clock":"time.perf_counter_ns","warmup":"one untimed 1x1 certificate",
            "timing_scope":{"construction":"complete make_certificate call including one-step baseline; excludes JSON output",
                            "checking":"complete verify_certificate call, including regeneration and JSON equality comparison",
                            "norm_only":"request normalization, I-alpha*A construction, one-step induced-norm product and exact budget search; shares these helpers with production, skips supports",
                            "direct_power":"complete direct_power_reference call including all smaller powers up to the first certificate or cap64"},
            "serialization_sizes":"compact JSON UTF-8 bytes, not peak RAM",
            "limitations":["one environment; three sequential observations; no solver runtime study",
                           "reference_limit is an unsuccessful capped run, not a successful competing certificate",
                           "exact-input models; no physical-model accuracy or floating-point solver validation"],
            "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"observations":commands}
    args.output.with_suffix(".json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
