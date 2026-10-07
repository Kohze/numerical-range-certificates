"""Convenience API for consuming certificate JSON."""
import json
from .core import CertificateError, verify_certificate

def check_json(text, *, allow_conditional=False):
    try:
        certificate = json.loads(text)
    except (TypeError, ValueError, UnicodeError, RecursionError) as exc:
        raise CertificateError("invalid certificate JSON") from exc
    return verify_certificate(certificate, allow_conditional=allow_conditional)

