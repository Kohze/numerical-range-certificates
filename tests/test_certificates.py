"""Mathematical invariant, provenance, and certificate tampering tests."""
import copy
import json
import random
import unittest
from fractions import Fraction

from nrcert import CertificateError, make_certificate, power_witness_leq, sqrt_upper, verify_certificate
from nrcert.core import (
    _evaluate, _matrix, _norm_squared, _supports, _vertices, rational,
)

def coupled(n):
    return [[
        "1" if i == j and i % 2 == 0 else
        "11/10" if i == j else
        "-4/5" if j == i + 1 and i % 2 == 0 else
        "-1/100" if j == i + 1 else "0"
        for j in range(n)
    ] for i in range(n)]

def vector_numerical_range(A, u, v):
    denominator = sum(a * a + b * b for a, b in zip(u, v))
    real = sum(
        A[i][j] * (u[i] * u[j] + v[i] * v[j])
        for i in range(len(A)) for j in range(len(A))
    ) / denominator
    imag = sum(
        A[i][j] * (u[i] * v[j] - v[i] * u[j])
        for i in range(len(A)) for j in range(len(A))
    ) / denominator
    return real, imag

class ExactCertificates(unittest.TestCase):
    def test_square_root_outward_and_minimal(self):
        generator = random.Random(17)
        for _ in range(200):
            x = Fraction(generator.randrange(100000), generator.randrange(1, 100000))
            upper = sqrt_upper(x, 8)
            self.assertGreaterEqual(upper * upper, x)
            if upper:
                self.assertLess((upper - Fraction(1, 256)) ** 2, x)

    def test_square_root_exact_and_zero(self):
        self.assertEqual(sqrt_upper(Fraction(9, 16)), Fraction(3, 4))
        self.assertEqual(sqrt_upper(0), 0)

    def test_float_rejected(self):
        with self.assertRaises(CertificateError):
            make_certificate({"matrix": [[1.1]]})

    def test_invalid_shapes_and_values(self):
        for A in ([], [[1, 2]], [[True]], [[None]], [[1], [2]]):
            with self.subTest(A=A), self.assertRaises(CertificateError):
                make_certificate({"matrix": A})

    def test_invalid_settings(self):
        for settings in (
            {"alpha": 0}, {"alpha": -1}, {"tolerance": 1},
            {"tolerance": 0}, {"sqrt_bits": True},
            {"direction_resolution": 0}, {"mode": "sharp"},
            {"unexpected": 1}, {"kind": "complex"},
        ):
            with self.subTest(settings=settings), self.assertRaises(CertificateError):
                make_certificate({"matrix": [[1]], **settings})

    def test_support_contains_exact_complex_rayleigh_quotients(self):
        generator = random.Random(812)
        for n in range(1, 5):
            for _ in range(12):
                A = [[Fraction(generator.randrange(-7, 8), 5)
                      for _ in range(n)] for _ in range(n)]
                supports = _supports(A, 2, 16)
                for _ in range(15):
                    u = [generator.randrange(-3, 4) for _ in range(n)]
                    v = [generator.randrange(-3, 4) for _ in range(n)]
                    if not any(u + v):
                        u[0] = 1
                    x, y = vector_numerical_range(A, u, v)
                    self.assertTrue(all(a*x + b*y <= h for a, b, h, _ in supports))

    def test_every_vertex_satisfies_every_halfplane(self):
        A = _matrix(coupled(8))
        supports = _supports(A, 4, 24)
        vertices = _vertices(supports)
        self.assertTrue(vertices)
        self.assertTrue(all(
            a*x + b*y <= h for x, y in vertices for a, b, h, _ in supports
        ))

    def test_point_polygon_and_zero_residual(self):
        certificate = make_certificate({"matrix": [[2]], "alpha": "1/2"})
        self.assertEqual(certificate["vertices"], [["2", "0"]])
        self.assertEqual(certificate["evidence"]["iterations"], 1)
        self.assertEqual(certificate["evidence"]["bound_squared"]["base"], "0")
        self.assertTrue(power_witness_leq(certificate["evidence"]["bound_squared"], 0))

    def test_segment_polygon(self):
        certificate = make_certificate({"matrix": [[1, 0], [0, 2]], "alpha": "2/3"})
        self.assertEqual(certificate["vertices"], [["1", "0"], ["2", "0"]])
        self.assertEqual(certificate["evidence"]["q_upper_squared"], "1/9")

    def test_coupled_matrix_beats_induced_norm_baseline(self):
        certificate = make_certificate({
            "matrix": coupled(8), "tolerance": "1/100000000",
        })
        evidence = certificate["evidence"]
        self.assertEqual(evidence["iterations"], 31)
        self.assertEqual(evidence["induced_norm_baseline"]["iterations"], 113)
        self.assertEqual(evidence["induced_norm_baseline"]["q_upper_squared"], "18/25")
        self.assertEqual(evidence["best_method"], "numerical_range")
        self.assertTrue(power_witness_leq(
            evidence["bound_squared"], evidence["tolerance_squared"],
        ))
        self.assertFalse(power_witness_leq(
            evidence["previous_bound_squared"], evidence["tolerance_squared"],
        ))

    def test_bidiagonal_control_favors_baseline(self):
        A = [["1" if i == j else "1/5" if j == i+1 else "0"
              for j in range(16)] for i in range(16)]
        evidence = make_certificate({"matrix": A, "tolerance": "1/100000000"})["evidence"]
        self.assertEqual(evidence["iterations"], 13)
        self.assertEqual(evidence["induced_norm_baseline"]["iterations"], 12)
        self.assertEqual(evidence["best_method"], "induced_norm")

    def test_geometric_failure_is_reported(self):
        certificate = make_certificate({"matrix": [[1, 2], [0, 1]]})
        self.assertEqual(certificate["evidence"]["status"], "no_geometric_contraction")
        self.assertIsNone(certificate["evidence"]["iterations"])
        self.assertFalse(verify_certificate(certificate)["target_certified"])

    def test_iteration_limit_is_reported(self):
        certificate = make_certificate({
            "matrix": [[1, 0], [0, 2]], "alpha": "2/3", "max_iterations": 1,
        })
        self.assertEqual(certificate["evidence"]["status"], "iteration_limit")
        self.assertIsNone(certificate["evidence"]["iterations"])

    def test_conditional_requires_explicit_opt_in(self):
        certificate = make_certificate({
            "kind": "polynomial", "matrix": [[1]], "coefficients": [1],
            "mode": "conditional_two",
        })
        with self.assertRaises(CertificateError):
            verify_certificate(certificate)
        result = verify_certificate(certificate, allow_conditional=True)
        self.assertTrue(result["conditional"])

    def test_deterministic_and_json_round_trip(self):
        request = {"matrix": coupled(8)}
        first = make_certificate(request)
        self.assertEqual(first, make_certificate(request))
        self.assertTrue(verify_certificate(json.loads(json.dumps(first)))["valid"])

    def test_tampered_evidence_rejected(self):
        certificate = make_certificate({"matrix": [[1, 0], [0, 2]], "alpha": "2/3"})
        for mutate in (
            lambda c: c["vertices"][0].__setitem__(0, "0"),
            lambda c: c.__setitem__("constant_upper", "1"),
            lambda c: c["support_evidence"][0].__setitem__("upper", "0"),
            lambda c: c["evidence"].__setitem__("iterations", 1),
            lambda c: c.__setitem__("algorithm", "trust-me"),
            lambda c: c.__setitem__("verified", True),
        ):
            changed = copy.deepcopy(certificate)
            mutate(changed)
            with self.assertRaises(CertificateError):
                verify_certificate(changed)

    def test_boolean_and_float_tampering_rejected(self):
        certificate = make_certificate({"matrix": [[1]]})
        for replacement in (False, 0.0):
            changed = copy.deepcopy(certificate)
            normal = changed["support_evidence"][1]["normal"]
            zero_index = normal.index(0)
            normal[zero_index] = replacement
            with self.assertRaises(CertificateError):
                verify_certificate(changed)

    def test_constant_and_zero_polynomials(self):
        for coefficients, expected in (([0], "0"), (["3/2"], "3/2")):
            evidence = make_certificate({
                "kind": "polynomial", "matrix": [[1]], "coefficients": coefficients,
            })["evidence"]
            self.assertEqual(evidence["scalar_supremum_upper"], expected)
            self.assertEqual(evidence["derivative_upper"], "0")

    def test_polynomial_segment_mesh(self):
        coarse = make_certificate({
            "kind": "polynomial", "matrix": [[1,0], [0,2]],
            "coefficients": [1, "-4/3", "4/9"], "samples_per_edge": 16,
        })
        fine = make_certificate({
            **coarse["request"], "samples_per_edge": 32,
        })
        self.assertEqual(fine["evidence"]["sample_count"], 33)
        self.assertEqual(fine["evidence"]["max_sample_squared"], "1/81")
        self.assertLess(
            rational(fine["evidence"]["scalar_supremum_upper"]),
            rational(coarse["evidence"]["scalar_supremum_upper"]),
        )

    def test_polynomial_bound_covers_rational_interior_points(self):
        request = {
            "kind": "polynomial", "matrix": coupled(4),
            "coefficients": [1, -2, 1], "samples_per_edge": 32,
        }
        certificate = make_certificate(request)
        vertices = [[rational(c) for c in z] for z in certificate["vertices"]]
        upper = rational(certificate["evidence"]["scalar_supremum_upper"])
        generator = random.Random(29)
        for _ in range(100):
            weights = [generator.randrange(1, 100) for _ in vertices]
            total = sum(weights)
            z = tuple(sum(w*v[k] for w,v in zip(weights,vertices))/total for k in (0,1))
            self.assertLessEqual(_norm_squared(_evaluate([Fraction(1),-2,1],z)),upper**2)

    def test_slow_contraction_uses_compact_power_witness(self):
        certificate = make_certificate({
            "matrix": [["1/100"]], "tolerance": "1/100000000",
        })
        self.assertEqual(certificate["evidence"]["iterations"], 1902)
        self.assertEqual(certificate["evidence"]["induced_norm_baseline"]["iterations"], 1833)
        encoded = json.dumps(certificate)
        self.assertLess(len(encoded), 5000)
        self.assertTrue(verify_certificate(json.loads(encoded))["target_certified"])
        self.assertTrue(power_witness_leq(
            certificate["evidence"]["bound_squared"], "1/10000000000000000",
        ))

    def test_invalid_power_witness_and_output_budget(self):
        for witness in (
            {}, {"coefficient":"-1", "base":"1", "exponent":1},
            {"coefficient":"1", "base":"1", "exponent":True},
        ):
            with self.assertRaises(CertificateError):
                power_witness_leq(witness, 1)
        with self.assertRaises(CertificateError):
            make_certificate({
                "kind":"polynomial", "matrix":[["3/2"]],
                "coefficients":["1/" + str(1 << 8000), "1/" + str((1 << 8000) + 1)],
            })

    def test_polynomial_bad_coefficients(self):
        for coefficients in ([], [1.1], "1,2"):
            with self.assertRaises(CertificateError):
                make_certificate({"kind":"polynomial", "matrix":[[1]],
                                  "coefficients":coefficients})

if __name__ == "__main__":
    unittest.main()

