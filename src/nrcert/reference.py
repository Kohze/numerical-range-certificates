"""Exact direct-power reference bounds for Richardson's residual.

This reference constructs matrix powers and applies
||B**m||_2**2 <= ||B**m||_1 ||B**m||_infinity.
It is intentionally separate from the numerical-range certificate producer.
The linear search checks every smaller positive power; no monotonicity of
matrix-power norms is assumed. A search-limit result does not prove divergence.
"""
from fractions import Fraction
from nrcert.core import CertificateError, _normalise, _str, rational

def _multiply(left, right):
    """Generic exact matrix multiplication, not a bidiagonal recurrence."""
    columns = list(zip(*right))
    return [
        [
            sum((x * y for x, y in zip(row, column) if x and y), Fraction(0))
            for column in columns
        ]
        for row in left
    ]

def _norm_product(A):
    row_norm = max(sum(abs(x) for x in row) for row in A)
    column_norm = max(sum(abs(x) for x in column) for column in zip(*A))
    return row_norm * column_norm

def direct_power_reference(request, max_power=64):
    """Return the first Richardson power passing the norm-product criterion.

    Accept the same validated Richardson fields as the certificate producer.
    Geometry and theorem-mode fields do not alter this constant-one reference.
    max_power is its separate finite scan limit, not request.max_iterations.
    """
    if not isinstance(request, dict):
        raise CertificateError("reference request must be an object")
    if type(max_power) is not int or not 1 <= max_power <= 128:
        raise CertificateError("reference max_power must lie in [1,128]")
    if request.get("kind", "richardson") != "richardson":
        raise CertificateError("direct-power reference supports Richardson requests only")
    request, A = _normalise(request)
    alpha = rational(request["alpha"], "alpha")
    tolerance = rational(request["tolerance"], "tolerance")
    n = len(A)
    B = [[Fraction(int(i == j)) - alpha * A[i][j] for j in range(n)] for i in range(n)]
    power = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    checked = []
    previous = Fraction(1)
    for m in range(1, max_power + 1):
        power = _multiply(power, B)
        product = _norm_product(power)
        checked.append(_str(product))
        if product <= tolerance ** 2:
            return {
                "request": request,
                "method": "exact_direct_power_induced_norm",
                "status": "target_certified",
                "first_certified_power": m,
                "bound_squared": _str(product),
                "previous_bound_squared": _str(previous),
                "checked_bound_squared": checked,
                "max_power": max_power,
            }
        previous = product
    return {
        "request": request,
        "method": "exact_direct_power_induced_norm",
        "status": "reference_limit",
        "first_certified_power": None,
        "bound_squared": None,
        "previous_bound_squared": None,
        "checked_bound_squared": checked,
        "max_power": max_power,
    }


