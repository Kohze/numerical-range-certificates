"""Exact numerical-range certificates for polynomial matrix calculations."""
from .core import (
    CertificateError,
    make_certificate,
    power_witness_leq,
    sqrt_upper,
    verify_certificate,
)

__version__ = "0.1.0"
__all__ = [
    "CertificateError", "make_certificate", "power_witness_leq",
    "sqrt_upper", "verify_certificate",
]
