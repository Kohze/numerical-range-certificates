"""Adversarial regressions from the first internal software referee round."""
import contextlib
import io
import json
import unittest
from fractions import Fraction
from unittest.mock import patch

from nrcert import (
    CertificateError, make_certificate, power_witness_leq, verify_certificate,
)
from nrcert.checker import check_json
from nrcert.cli import main
from nrcert.core import rational
from nrcert.reference import direct_power_reference


def invoke(text, command="certify"):
    output, errors = io.StringIO(), io.StringIO()
    with patch("sys.stdin", io.StringIO(text)):
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            code = main([command, "-"])
    return code, output.getvalue(), errors.getvalue()


class RefereeRegressions(unittest.TestCase):
    def test_generated_large_factors_are_comparable(self):
        denominator = 2**4500
        certificate = make_certificate({
            "matrix": [[Fraction(denominator+1, 2*denominator)]],
            "tolerance": "1/1000", "max_iterations": 64,
        })
        evidence = certificate["evidence"]
        self.assertGreater(len(evidence["bound_squared"]["base"]), 4096)
        self.assertGreater(Fraction(evidence["bound_squared"]["base"]).denominator.bit_length(), 8192)
        self.assertTrue(verify_certificate(certificate)["target_certified"])
        self.assertTrue(power_witness_leq(
            evidence["bound_squared"], evidence["tolerance_squared"],
        ))
        self.assertFalse(power_witness_leq(
            evidence["previous_bound_squared"], evidence["tolerance_squared"],
        ))

    def test_large_normalized_inputs_remain_recomputable(self):
        denominator = 2**7000
        for value in (
            Fraction(denominator+1, 2*denominator),
            "5"+"1"*2099+"e-2100",
        ):
            with self.subTest(value_type=type(value).__name__):
                certificate = make_certificate({
                    "kind": "polynomial", "matrix": [[value]], "coefficients": [0],
                })
                self.assertGreater(len(certificate["request"]["matrix"][0][0]), 4096)
                result = check_json(json.dumps(certificate))
                self.assertTrue(result["valid"])
                self.assertTrue(result["bound_certified"])

    def test_derived_factor_and_threshold_limits(self):
        witness = {"coefficient": "1", "base": "1", "exponent": 1}
        for field in ("coefficient", "base"):
            with self.subTest(field=field), self.assertRaises(CertificateError):
                power_witness_leq({**witness, field: 1 << 12000}, 1)
        with self.assertRaises(CertificateError):
            power_witness_leq(witness, 1 << 12000)

    def test_compact_exponents_rejected_before_fraction_allocation(self):
        arabic = "1e" + "".join(chr(0x660+int(ch)) for ch in "1000000000")
        fullwidth = "1e-" + "".join(chr(0xFF10+int(ch)) for ch in "1000000000")
        for text in ("1e1000000000", "1e-1000000000",
                     "0e1000000000", "1e+1_000_000_000", arabic, fullwidth):
            with self.subTest(text=text), patch("nrcert.core.Fraction", wraps=Fraction) as constructor:
                with self.assertRaisesRegex(CertificateError, "exponent"):
                    rational(text)
                constructor.assert_not_called()
        for text in ("1e1000000000", arabic, fullwidth):
            with self.subTest(derived=text), patch("nrcert.core.Fraction", wraps=Fraction) as constructor:
                with self.assertRaisesRegex(CertificateError, "exponent"):
                    power_witness_leq(
                        {"coefficient": "1", "base": text, "exponent": 1}, 1,
                    )
                # Coefficient 1 is parsed first. The enormous base is never parsed.
                self.assertEqual(constructor.call_count, 1)

    def test_ordinary_exact_scientific_decimals_remain_supported(self):
        for text, expected in (
            ("1e2", Fraction(100)), ("-2.5e-3", Fraction(-1, 400)),
            ("1e+03", Fraction(1000)), (".125", Fraction(1, 8)),
            ("123/456", Fraction(41, 152)),
            ("1e" + chr(0x661) + chr(0x660), Fraction(10**10)),
        ):
            with self.subTest(text=text):
                self.assertEqual(rational(text), expected)

    def test_unsupported_json_is_a_controlled_cli_failure(self):
        for text in ('{"matrix":[[' + "1"*5000 + "]]}",
                     "["*1200+"0"+"]"*1200):
            with self.subTest(length=len(text)):
                code, output, error = invoke(text)
                self.assertEqual(code, 1)
                self.assertEqual(output, "")
                self.assertIn("nrcert: error:", error)
                self.assertNotIn("Traceback", error)

    def test_unsupported_json_is_a_certificate_error(self):
        for text in ('{"matrix":[[' + "1"*5000 + "]]}",
                     "["*1200+"0"+"]"*1200, b"\xff"):
            with self.subTest(type=type(text).__name__):
                with self.assertRaises(CertificateError):
                    check_json(text)

    def test_richardson_default_has_unconditional_power_provenance(self):
        certificate = make_certificate({
            "matrix": [["1/100"]], "tolerance": "1/100000000",
        })
        self.assertEqual(certificate["constant_upper"], "2")
        self.assertIn("power inequality", certificate["assumption"])
        self.assertEqual(certificate["evidence"]["iterations"], 1902)
        self.assertFalse(verify_certificate(certificate)["conditional"])
        polynomial = make_certificate({
            "kind": "polynomial", "matrix": [[1]], "coefficients": [1],
        })
        self.assertGreater(Fraction(polynomial["constant_upper"]), 2)
        self.assertIn("Crouzeix", polynomial["assumption"])

    def test_explicit_crouzeix_palencia_power_comparator(self):
        certificate = make_certificate({
            "matrix": [["1/100"]], "tolerance": "1/100000000",
            "mode": "crouzeix_palencia",
        })
        self.assertGreater(Fraction(certificate["constant_upper"]), 2)
        self.assertEqual(certificate["evidence"]["iterations"], 1921)
        self.assertFalse(verify_certificate(certificate)["conditional"])

    def test_conditional_two_power_request_has_clear_migration_error(self):
        request = {"matrix": [[1]], "mode": "conditional_two"}
        with self.assertRaisesRegex(CertificateError, "use mode 'established'"):
            make_certificate(request)
        code, output, error = invoke(json.dumps(request))
        self.assertEqual(code, 1)
        self.assertEqual(output, "")
        self.assertIn("unconditional constant 2", error)

    def test_reference_strict_request_and_recorded_normalized_parameters(self):
        for request in (
            {"kind": "polynomial", "matrix": [[1]], "coefficients": [7]},
            {"matrix": [["1/2"]], "tolerence": "1/100000000"},
        ):
            with self.subTest(request=request), self.assertRaises(CertificateError):
                direct_power_reference(request)
        request = {
            "matrix": [["1/2"]], "tolerance": "1/100000000",
            "mode": "crouzeix_palencia", "direction_resolution": 4,
        }
        result = direct_power_reference(request)
        self.assertEqual(result["first_certified_power"], 27)
        self.assertEqual(result["request"]["tolerance"], request["tolerance"])
        self.assertEqual(result["request"]["alpha"], "1")
        self.assertEqual(result["request"]["kind"], "richardson")


if __name__ == "__main__":
    unittest.main()

