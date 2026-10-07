# Numerical Range Certificates

`nrcert` produces exact-rational upper bounds for polynomial matrix computations and iteration budgets for Richardson's method. It is intended for researchers who need a reviewable bound for a small or moderate real rational matrix, together with enough evidence to recompute it.

The default `established` mode uses the classical unconditional numerical-radius power bound `C = 2` for Richardson iteration and the Crouzeix--Palencia bound `C = 1 + sqrt(2)` for general polynomials. `crouzeix_palencia` explicitly selects the latter theorem for either kind. For general polynomials only, an opt-in `conditional_two` mode reports what the sharp-Crouzeix hypothesis `C = 2` would imply. The package does not prove that hypothesis or certify the OpenAI repository's proofs.

## Install and run

Python 3.10 or later is required. There are no runtime dependencies.

```sh
python -m pip install .
python -m unittest discover -s tests -v
nrcert certify examples/coupled32.json -o certificate.json
nrcert check certificate.json
```

Inputs accept integers and exact rational or decimal strings such as `"4/5"` or `"0.8"`. JSON floating-point numbers, including `0.8` without quotes, are rejected. A certificate applies to the supplied exact matrix.

An input for a normal control is:

```json
{
  "kind": "richardson",
  "matrix": [[1, 0], [0, 2]],
  "alpha": "2/3",
  "tolerance": "1/100000000",
  "mode": "established",
  "direction_resolution": 1,
  "sqrt_bits": 40,
  "max_iterations": 10000
}
```

For Richardson `x_(m+1) = x_m + alpha (b - A x_m)`, the exact residual satisfies `r_m = (I - alpha A)^m r_0`. The certificate bounds `||r_m||_2 / ||r_0||_2` for every initial residual. It does not certify the residual produced by a rounded implementation.

The normal example yields 18 steps from the universal numerical-range theorem and 17 from the induced-norm baseline. The certificate reports both and selects the smaller successful budget among these two emitted bounds. A separate exact direct-power reference can give a stronger budget. A normal matrix also has a spectral-theorem bound with constant one; this example is a control.

The nonnormal `coupled32.json` fixture has alternating diagonal entries `1, 11/10` and alternating superdiagonal entries `-4/5, -1/100`. At tolerance `10^-8`, eight directions certify 31 steps; the cheap one-step induced-norm bound gives 113, and the stronger exact direct-power reference gives 11. Both unconditional Richardson theorem routes give 31 on this fixture. Thirty-two directions improve the default budget to 29, which still loses to the direct-power reference. These are certified upper budgets, not observed iteration counts, industrial validation, or evidence of overall numerical-range method superiority.

## Polynomial norm bounds

Coefficients are listed from the constant term upward:

```sh
nrcert certify examples/polynomial_segment.json -o polynomial-certificate.json
nrcert check polynomial-certificate.json
```

That request uses `p(z) = (1 - 2z/3)^2` on `diag(1,2)`. The producer meshes polygon edges with exact rational points, evaluates the polynomial exactly, and adds a derivative-based mesh allowance. Points and line segments are handled directly. The result bounds `||p(A)||_2` in exact arithmetic. A coefficient/induced-norm baseline is also reported.

## Python API

```python
from nrcert import make_certificate, power_witness_leq, verify_certificate

certificate = make_certificate({
    "matrix": [[1, 0], [0, 2]],
    "alpha": "2/3",
    "tolerance": "1/100000000",
})
report = verify_certificate(certificate)
evidence = certificate["evidence"]
assert report["target_certified"]
assert power_witness_leq(
    evidence["bound_squared"], evidence["tolerance_squared"]
)
```

Squared Richardson bounds are stored as compact power witnesses:
`{"coefficient": c_squared, "base": q_squared, "exponent": m}`.
Their value is `c_squared * q_squared**m`. This avoids serializing enormous expanded integer powers when contraction is slow. Every iteration comparison still uses exact integer/rational arithmetic.

The checker recomputes supports, vertices, samples, bounds, and the first successful budget for each of the two certificate bounds from the normalized request. `comparison_methods` and `comparison_scope` explicitly restrict `best_certified_iterations` to those two bounds; the method label `induced_norm` denotes the one-step bound. This is not a global optimality claim. It rejects evidence inconsistent with the declared request and unknown evidence. The caller must compare the declared matrix with the intended input; the checker does not authenticate where an input came from. Producer and checker share the same implementation; recomputation is not an independent implementation audit or a formal proof of the software.

## Mathematical basis

For integer normal `(a,b)`, write `S = (A + A.T)/2` and `K = (A - A.T)/2`. The Hermitian support matrix is `a S - i b K`. Gershgorin's theorem gives an upper support using row diagonals plus off-diagonal magnitudes. Each magnitude is rounded outward to a multiple of `2^-sqrt_bits` using integer square roots. Intersections of these valid halfplanes give an outer polygon for the complex numerical range, even when A is real.

For `p(z) = (1-alpha z)^m`, the largest `|1-alpha z|` over a compact polygon is attained at a vertex. The classical inequalities `w(B^m) <= w(B)^m` and `||B||_2 <= 2 w(B)` give `||p(A)||_2 <= 2 q^m` for every positive integer m, without a sharp-Crouzeix assumption. See [Tadmor (2025), section 3.2](https://www.math.umd.edu/~tadmor/pub/linear-stability/Tadmor%20RK%20are%20stable%20CPAM2025.pdf). For other polynomials, the boundary mesh adds `delta * sum_k k |c_k| rho^(k-1)`, with outward rational radius and mesh bounds.

Support-line numerical-range construction is established work; see [Johnson (1978)](https://doi.org/10.1137/0715039) and [Chebfun's field-of-values example](https://www.chebfun.org/examples/linalg/FieldOfValues.html). This package's focus is exact rational evidence, explicit theorem provenance, and recomputation. The default general-polynomial constant comes from [Crouzeix and Palencia (2017)](https://doi.org/10.1137/17M1116672). The broader distinction between a numerical calculation and a rigorous enclosure is discussed by [Rump (2010)](https://doi.org/10.1017/S096249291000005X).

## Assumptions and failure modes

- For general polynomials, `conditional_two` must be requested explicitly. Both CLI construction and checking additionally require `--allow-conditional`. A valid conditional certificate remains conditional. Richardson requests reject that mode because their constant-two route is unconditional.
- `no_geometric_contraction` means this geometric bound has `q >= 1`. It does not prove Richardson diverges. The shifted nilpotent control has this status although its exact residual vanishes after two steps.
- `iteration_limit` means no budget was certified within the requested search limit.
- Refining directions reduces angular discretization but does not remove Gershgorin slack. There is no unconditional convergence claim to the true numerical range.
- The present implementation accepts real rational matrices only, dimension at most 64, direction resolution 1--8 (8--64 directions), polynomial degree at most 128, at most 128 mesh segments per edge, and at most 100000 iterations. Exact arithmetic can be expensive.
- Input rational strings have at most 8192 characters and at most 8192-bit numerators/denominators. Derived rational factors have an 8192-character and 12000-bit budget. Decimal exponents are checked before conversion, including supported Unicode digits. Expanded output rationals have a 12000-bit budget; excessive precision, degree, or coefficient growth raises `CertificateError`. Factored iteration bounds avoid serializing expanded powers. These size caps do not guarantee low computation cost or make the API a hardened hostile-input service.
- These are exact-model certificates. Input uncertainty, matrix assembly error, and solver roundoff require additional bounds. A float converted to a string is not a validation of its correspondence to the intended model.

Exit codes are 0 for a produced/recomputed bound, 1 for invalid input or evidence, and 2 for valid geometric evidence without a successful bound. JSON is written to standard output unless `--output` is supplied; `-` reads standard input.

## Paper and reproducibility

The companion [numerical paper](paper/numerical-range-certificates.pdf) presents the finite-input theorem, error decomposition, analytic controls and comparative studies. Its [LaTeX source](paper/main.tex), [bibliography](paper/references.bib) and [reference-verification ledger](audit/citation-verification.json) are included. The theorem permits matrix-valued coefficients; the shipped prototype implements real-rational scalar polynomials.

This alpha prototype provides reviewable evidence and quantified diagnostics. It demonstrates no general solver-performance advantage: on the featured coupled matrix, the exact-power reference is both stronger and cheaper than numerical-range construction in the recorded observations. Producer/checker recomputation, independently implemented finite-study checks, and a formal verification of all software are distinct levels of assurance.

After installing the package, reproduce or verify the public data with:

```sh
python -B experiments/reproduce.py
python -B experiments/direct_power_baseline.py examples/coupled32.json
python -B reproduce/check_saved_results.py
python -B experiments/independent_ellipse_audit.py
```

`check_saved_results.py` compares all 30 comparison rows and 144 exact ellipse-study rows, checks saved certificates, and verifies the recorded timing medians. It does not require NumPy and does not overwrite the retained results. To generate a new study:

```sh
python -B experiments/ellipse_error_budget.py --output output/ellipse-error-budget.csv
python -B experiments/certificate_cost_study.py --output output/certificate-cost-study.csv
```

The cost JSON retains all original observations and their environment/source provenance. New timings are descriptive measurements on the machine executing the script and need not match the paper. Optional `python experiments/exploratory_norms.py` requires NumPy; those floating-point spectral norms are diagnostics.

For a separately installed wheel, `python -B reproduce/verify_package.py` executes the tests, comparison reproduction and six CLI examples, and writes portable reports under `audit/`. Build the paper with `python -B reproduce/build_paper.py` (requires `pdflatex` and `bibtex`); the rebuilt PDF goes to `output/pdf/`.

## Project and citation

Project: [Kohze/numerical-range-certificates](https://github.com/Kohze/numerical-range-certificates). Author and maintainer: Robin Gounder, Vaionex Corporation; [robin@gounder.com](mailto:robin@gounder.com). See [CITATION.cff](CITATION.cff) and [CONTRIBUTING.md](CONTRIBUTING.md).

Version 0.1.0 is an alpha release, available with data and installation files from the [tagged release](https://github.com/Kohze/numerical-range-certificates/releases/tag/v0.1.0). The recorded preparation environment was Python 3.11.5. CI checks Python 3.10–3.13; completed results appear under [Actions](https://github.com/Kohze/numerical-range-certificates/actions).

## Licensing

The software, tests, experiment/reproduction programs and package documentation are MIT-licensed: see [LICENSE](LICENSE). The manuscript source and PDF in `paper/` retain their author's copyright under [paper/LICENSE](paper/LICENSE); the software licence does not grant a manuscript reuse licence.
