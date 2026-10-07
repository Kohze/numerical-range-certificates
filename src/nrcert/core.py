"""Exact-rational numerical-range certificates for real rational matrices.

No floating-point arithmetic is used in certificate construction or checking.
Conditional constant-two bounds are never selected implicitly.
"""
from fractions import Fraction
from math import isqrt
import json
import re

SCHEMA = "nrcert/1"
ALGORITHM = "rational-gershgorin-polygon-v1"
MAX_DIMENSION = 64
MAX_RESOLUTION = 8
MAX_DEGREE = 128
MAX_SAMPLES_PER_EDGE = 128
MAX_ITERATIONS = 100000
MAX_EXPANDED_BITS = 12000
MAX_INPUT_CHARS = 8192
MAX_INPUT_BITS = 8192
MAX_DERIVED_CHARS = 8192  # covers two 12000-bit decimal integers and a slash
_DECIMAL_EXPONENT = re.compile(r"[eE]([+-]?\d(?:_?\d)*)$")
POWER_ASSUMPTION = (
    "Numerical-radius power inequality: C = 2 for Richardson powers "
    "(m >= 1), unconditional"
)
ESTABLISHED_ASSUMPTION = "Crouzeix--Palencia (2017): C = 1 + sqrt(2)"
CONDITIONAL_ASSUMPTION = (
    "Hypothesis H2: ||p(A)||_2 <= 2 sup_{z in W(A)} |p(z)| "
    "for every polynomial p; this checker does not establish H2"
)

class CertificateError(ValueError):
    """Invalid request, evidence, or unsupported mathematical assumptions."""

def _integer(value, name, low, high):
    if type(value) is not int or not low <= value <= high:
        raise CertificateError(f"{name} must be an integer in [{low}, {high}]")
    return value

def _bounded_rational(value, name, max_chars, max_bits):
    """Apply lexical/exponent limits before Fraction can allocate large powers."""
    if type(value) is int:
        if value.bit_length() > max_bits:
            raise CertificateError(f"{name} exceeds the supported rational size")
        result = Fraction(value)
    elif isinstance(value, str) and len(value) <= max_chars:
        # Fraction accepts scientific decimals. A short exponent can otherwise
        # allocate arbitrarily large integers before the final bit-size check.
        match = _DECIMAL_EXPONENT.search(value.strip())
        if match:
            # Fraction also accepts Unicode decimal digits. Normalize one
            # bounded digit at a time so no exponent spelling bypasses this guard.
            digits = "".join(
                str(int(ch)) for ch in match.group(1).lstrip("+-") if ch != "_"
            ).lstrip("0") or "0"
            limit = str(max_bits)
            if len(digits) > len(limit) or (len(digits) == len(limit) and digits > limit):
                raise CertificateError(f"{name} decimal exponent exceeds the parsing budget")
        try:
            result = Fraction(value)
        except (ValueError, ZeroDivisionError, OverflowError) as exc:
            raise CertificateError(f"{name} is not an exact rational") from exc
    elif isinstance(value, Fraction):
        result = value
    else:
        raise CertificateError(f"{name} must be an integer or exact string")
    if max(result.numerator.bit_length(), result.denominator.bit_length()) > max_bits:
        raise CertificateError(f"{name} exceeds the supported rational size")
    return result

def rational(value, name="rational"):
    """Read bounded exact inputs; never accept floats."""
    return _bounded_rational(value, name, MAX_INPUT_CHARS, MAX_INPUT_BITS)

def _derived_rational(value, name):
    """Read factors produced within the larger expanded-output budget."""
    return _bounded_rational(value, name, MAX_DERIVED_CHARS, MAX_EXPANDED_BITS)

def _str(x):
    if max(x.numerator.bit_length(), x.denominator.bit_length()) > MAX_EXPANDED_BITS:
        raise CertificateError(
            "expanded rational exceeds the 12000-bit output budget; "
            "reduce precision, coefficient size, or polynomial degree"
        )
    try:
        return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
    except ValueError as exc:
        raise CertificateError("Python's integer-to-string budget is too small") from exc

def _power_witness(coefficient, base, exponent):
    return {"coefficient": _str(coefficient), "base": _str(base), "exponent": exponent}

def power_witness_leq(witness, threshold):
    """Compare a factored squared bound to a rational threshold exactly."""
    if not isinstance(witness, dict) or set(witness) != {"coefficient", "base", "exponent"}:
        raise CertificateError("invalid power witness")
    coefficient = _derived_rational(witness["coefficient"], "witness coefficient")
    base = _derived_rational(witness["base"], "witness base")
    exponent = _integer(witness["exponent"], "witness exponent", 0, MAX_ITERATIONS)
    if coefficient < 0 or base < 0:
        raise CertificateError("power witness must have nonnegative factors")
    return coefficient * base ** exponent <= _derived_rational(threshold, "threshold")

def sqrt_upper(x, bits=40):
    """Smallest multiple of 2**(-bits) whose square is at least x."""
    x = x if isinstance(x, Fraction) else rational(x)
    bits = _integer(bits, "sqrt_bits", 8, 256)
    if x < 0:
        raise CertificateError("cannot enclose the square root of a negative value")
    scale = 1 << bits
    k = isqrt((x.numerator * scale * scale) // x.denominator)
    if k * k * x.denominator < x.numerator * scale * scale:
        k += 1
    return Fraction(k, scale)

def _matrix(value):
    if not isinstance(value, list) or not value:
        raise CertificateError("matrix must be a nonempty square list")
    n = len(value)
    if n > MAX_DIMENSION:
        raise CertificateError(f"matrix dimension exceeds {MAX_DIMENSION}")
    if any(not isinstance(row, list) or len(row) != n for row in value):
        raise CertificateError("matrix must be square")
    return [[rational(x, "matrix entry") for x in row] for row in value]

def _directions(resolution):
    k = resolution
    # Every direction on the boundary of the integer square, without duplicates.
    return sorted(
        {(k, j) for j in range(-k, k + 1)}
        | {(-k, j) for j in range(-k, k + 1)}
        | {(j, k) for j in range(-k + 1, k)}
        | {(j, -k) for j in range(-k + 1, k)}
    )

def _supports(A, resolution, bits):
    """Exact rational upper bounds for lambda_max(a*S - i*b*K)."""
    n = len(A)
    supports = []
    for a, b in _directions(resolution):
        row_upper = []
        for i in range(n):
            value = a * A[i][i]
            for j in range(n):
                if j != i:
                    s = (A[i][j] + A[j][i]) / 2
                    k = (A[i][j] - A[j][i]) / 2
                    value += sqrt_upper((a * s) ** 2 + (b * k) ** 2, bits)
            row_upper.append(value)
        supports.append((a, b, max(row_upper), row_upper))
    return supports

def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

def _hull(points):
    points = sorted(set(points))
    if len(points) <= 1:
        return points
    lower, upper = [], []
    for p in points:
        while len(lower) >= 2 and _cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(points):
        while len(upper) >= 2 and _cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def _vertices(supports):
    points = []
    for i, (a, b, h, _) in enumerate(supports):
        for c, d, g, _ in supports[i + 1:]:
            determinant = a * d - b * c
            if determinant:
                x = (h * d - b * g) / determinant
                y = (a * g - h * c) / determinant
                if all(u * x + v * y <= w for u, v, w, _ in supports):
                    points.append((x, y))
    hull = _hull(points)
    if not hull:
        raise CertificateError("support polygon has no vertices")
    return hull

def _edges(vertices):
    if len(vertices) == 1:
        return [(vertices[0], vertices[0])]
    if len(vertices) == 2:
        return [(vertices[0], vertices[1])]
    return list(zip(vertices, vertices[1:] + vertices[:1]))

def _complex_mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])

def _evaluate(coefficients, z):
    result = (Fraction(0), Fraction(0))
    for coefficient in reversed(coefficients):
        result = _complex_mul(result, z)
        result = (result[0] + coefficient, result[1])
    return result

def _norm_squared(z):
    return z[0] ** 2 + z[1] ** 2

def _constant(mode, bits, kind):
    if mode == "established" and kind == "richardson":
        return Fraction(2), POWER_ASSUMPTION
    if mode in ("established", "crouzeix_palencia"):
        return 1 + sqrt_upper(Fraction(2), bits), ESTABLISHED_ASSUMPTION
    if mode == "conditional_two" and kind == "polynomial":
        return Fraction(2), CONDITIONAL_ASSUMPTION
    if mode == "conditional_two" and kind == "richardson":
        raise CertificateError(
            "Richardson powers have unconditional constant 2; use mode 'established'"
        )
    raise CertificateError(
        "mode must be 'established', 'crouzeix_palencia', or "
        "'conditional_two' (polynomial only)"
    )

def _normalise(request):
    if not isinstance(request, dict):
        raise CertificateError("request must be a JSON object")
    if any(not isinstance(key, str) for key in request):
        raise CertificateError("request keys must be strings")
    kind = request.get("kind", "richardson")
    common = {"kind", "matrix", "mode", "sqrt_bits", "direction_resolution"}
    allowed = common | (
        {"alpha", "tolerance", "max_iterations"} if kind == "richardson"
        else {"coefficients", "samples_per_edge"} if kind == "polynomial" else set()
    )
    if kind not in ("richardson", "polynomial"):
        raise CertificateError("kind must be 'richardson' or 'polynomial'")
    extra = set(request) - allowed
    if extra:
        raise CertificateError(f"unknown request fields: {', '.join(sorted(extra))}")
    A = _matrix(request.get("matrix"))
    bits = _integer(request.get("sqrt_bits", 40), "sqrt_bits", 8, 256)
    resolution = _integer(
        request.get("direction_resolution", 1), "direction_resolution",
        1, MAX_RESOLUTION,
    )
    mode = request.get("mode", "established")
    _constant(mode, bits, kind)
    result = {
        "kind": kind,
        "matrix": [[_str(x) for x in row] for row in A],
        "mode": mode,
        "sqrt_bits": bits,
        "direction_resolution": resolution,
    }
    if kind == "richardson":
        alpha = rational(request.get("alpha", "1"), "alpha")
        tolerance = rational(request.get("tolerance", "1/1000000"), "tolerance")
        if alpha <= 0:
            raise CertificateError("alpha must be positive")
        if not 0 < tolerance < 1:
            raise CertificateError("tolerance must lie strictly between zero and one")
        result.update({
            "alpha": _str(alpha),
            "tolerance": _str(tolerance),
            "max_iterations": _integer(
                request.get("max_iterations", 10000), "max_iterations",
                1, MAX_ITERATIONS,
            ),
        })
    else:
        raw = request.get("coefficients")
        if not isinstance(raw, list) or not 1 <= len(raw) <= MAX_DEGREE + 1:
            raise CertificateError(f"coefficients must contain 1 to {MAX_DEGREE + 1} entries")
        coefficients = [rational(x, "coefficient") for x in raw]
        while len(coefficients) > 1 and coefficients[-1] == 0:
            coefficients.pop()
        result.update({
            "coefficients": [_str(x) for x in coefficients],
            "samples_per_edge": _integer(
                request.get("samples_per_edge", 16), "samples_per_edge",
                1, MAX_SAMPLES_PER_EDGE,
            ),
        })
    return result, A

def _induced_product(A):
    one = max(sum(abs(row[j]) for row in A) for j in range(len(A)))
    infinity = max(sum(abs(x) for x in row) for row in A)
    return one * infinity

def _minimum_iterations(q2, c2, tolerance2, maximum):
    """Find the first certified m using exact comparisons and binary search."""
    if q2 >= 1:
        return None, "no_geometric_contraction"
    if q2 == 0:
        return 1, "target_certified"
    def meets(m):
        return c2 * q2 ** m <= tolerance2
    lower, hi = 0, 1
    while hi < maximum and not meets(hi):
        lower, hi = hi, min(2 * hi, maximum)
    if not meets(hi):
        return None, "iteration_limit"
    lo = lower + 1
    while lo < hi:
        mid = (lo + hi) // 2
        if meets(mid):
            hi = mid
        else:
            lo = mid + 1
    return lo, "target_certified"

def _richardson_evidence(request, vertices, constant, A):
    alpha = rational(request["alpha"])
    tolerance = rational(request["tolerance"])
    candidates = [
        ((1 - alpha * x) ** 2 + (alpha * y) ** 2, (x, y))
        for x, y in vertices
    ]
    q2, maximizer = max(candidates)
    c2 = constant ** 2
    m, status = _minimum_iterations(q2, c2, tolerance ** 2, request["max_iterations"])
    residual = [
        [(Fraction(int(i == j)) - alpha * x) for j, x in enumerate(row)]
        for i, row in enumerate(A)
    ]
    baseline2 = _induced_product(residual)
    baseline_m, baseline_status = _minimum_iterations(
        baseline2, Fraction(1), tolerance ** 2, request["max_iterations"],
    )
    choices = [(count, method) for count, method in (
        (m, "numerical_range"), (baseline_m, "induced_norm"),
    ) if count is not None]
    best = min(choices) if choices else (None, None)
    evidence = {
        "status": status,
        "q_upper_squared": _str(q2),
        "q_upper": _str(sqrt_upper(q2, request["sqrt_bits"])),
        "maximizing_vertex": [_str(x) for x in maximizer],
        "iterations": m,
        "bound_squared": None if m is None else _power_witness(c2, q2, m),
        "previous_bound_squared": (
            None if m is None or m == 1 else _power_witness(c2, q2, m - 1)
        ),
        "tolerance_squared": _str(tolerance ** 2),
        "induced_norm_baseline": {
            "q_upper_squared": _str(baseline2),
            "iterations": baseline_m,
            "status": baseline_status,
            "bound_squared": None if baseline_m is None else _power_witness(Fraction(1), baseline2, baseline_m),
        },
        "comparison_methods": ["numerical_range", "induced_norm"],
        "comparison_scope": "Only numerical-range and one-step induced-norm bounds; excludes direct powers",
        "best_certified_iterations": best[0],
        "best_method": best[1],
    }
    return evidence

def _polynomial_evidence(request, vertices, constant, A):
    coefficients = [rational(x) for x in request["coefficients"]]
    segments = request["samples_per_edge"]
    radius2 = max(_norm_squared(z) for z in vertices)
    radius = sqrt_upper(radius2, request["sqrt_bits"])
    derivative = sum(
        k * abs(c) * radius ** (k - 1)
        for k, c in enumerate(coefficients) if k
    )
    maximum2, maximizer, sample_count, edge2 = Fraction(-1), None, 0, Fraction(0)
    for first, last in _edges(vertices):
        dx, dy = last[0] - first[0], last[1] - first[1]
        edge2 = max(edge2, dx ** 2 + dy ** 2)
        for j in range(segments + 1):
            z = (first[0] + dx * j / segments, first[1] + dy * j / segments)
            value2 = _norm_squared(_evaluate(coefficients, z))
            if value2 > maximum2:
                maximum2, maximizer = value2, z
            sample_count += 1
    sample_upper = sqrt_upper(maximum2, request["sqrt_bits"])
    delta = sqrt_upper(edge2, request["sqrt_bits"]) / (2 * segments)
    upper = sample_upper + delta * derivative
    matrix_upper = sqrt_upper(_induced_product(A), request["sqrt_bits"])
    baseline = sum(abs(c) * matrix_upper ** k for k, c in enumerate(coefficients))
    return {
        "status": "bound_certified",
        "sample_count": sample_count,
        "max_sample_squared": _str(maximum2),
        "maximizing_sample": [_str(x) for x in maximizer],
        "radius_upper": _str(radius),
        "derivative_upper": _str(derivative),
        "mesh_radius_upper": _str(delta),
        "scalar_supremum_upper": _str(upper),
        "matrix_norm_upper": _str(constant * upper),
        "matrix_norm_upper_squared": _str((constant * upper) ** 2),
        "induced_norm_coefficient_baseline": _str(baseline),
    }

def make_certificate(request):
    """Produce a deterministic, exact certificate from a request dictionary."""
    request, A = _normalise(request)
    bits = request["sqrt_bits"]
    constant, assumption = _constant(request["mode"], bits, request["kind"])
    supports = _supports(A, request["direction_resolution"], bits)
    vertices = _vertices(supports)
    evidence = (
        _richardson_evidence(request, vertices, constant, A)
        if request["kind"] == "richardson"
        else _polynomial_evidence(request, vertices, constant, A)
    )
    return {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "request": request,
        "assumption": assumption,
        "constant_upper": _str(constant),
        "support_evidence": [
            {"normal": [a, b], "upper": _str(h),
             "gershgorin_row_uppers": [_str(x) for x in rows]}
            for a, b, h, rows in supports
        ],
        "vertices": [[_str(x), _str(y)] for x, y in vertices],
        "evidence": evidence,
    }

def verify_certificate(certificate, *, allow_conditional=False):
    """Recompute every mathematical field; a stored 'verified' flag is ignored.

    This verifies the rational calculation, not the underlying Crouzeix theorem,
    the input matrix's correspondence to a physical system, or a rounded solver.
    """
    if not isinstance(certificate, dict):
        raise CertificateError("certificate must be a JSON object")
    request = certificate.get("request")
    if not isinstance(request, dict):
        raise CertificateError("certificate has no request object")
    conditional = request.get("kind") == "polynomial" and request.get("mode") == "conditional_two"
    if conditional and not allow_conditional:
        raise CertificateError("conditional H2 certificate requires explicit opt-in")
    expected = make_certificate(request)
    try:
        actual_json = json.dumps(certificate, sort_keys=True, separators=(",", ":"), allow_nan=False)
        expected_json = json.dumps(expected, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise CertificateError("certificate must contain finite JSON values") from exc
    if actual_json != expected_json:
        raise CertificateError("certificate does not match recomputed exact evidence")
    status = expected["evidence"]["status"]
    target = (
        expected["evidence"]["best_certified_iterations"] is not None
        if request["kind"] == "richardson" else False
    )
    return {
        "valid": True,
        "conditional": conditional,
        "status": status,
        "bound_certified": target or status == "bound_certified",
        "target_certified": target,
        "scope": (
            "exact-arithmetic ||r_m||_2 / ||r_0||_2 bound"
            if request["kind"] == "richardson" else "exact-arithmetic ||p(A)||_2 bound"
        ),
    }


