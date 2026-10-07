"""Regression checks for a stronger exact-power comparison."""
import unittest
from fractions import Fraction
from nrcert import CertificateError
from nrcert.reference import _multiply, direct_power_reference

def coupled(n):
    return [[
        "1" if i == j and i % 2 == 0 else "11/10" if i == j else
        "-4/5" if j == i+1 and i % 2 == 0 else "-1/100" if j == i+1 else "0"
        for j in range(n)
    ] for i in range(n)]

class DirectPowerReference(unittest.TestCase):
    def test_generic_multiplication(self):
        matrix = [[Fraction(1), Fraction(2)], [Fraction(3), Fraction(4)]]
        self.assertEqual(_multiply(matrix,matrix), [[7,10],[15,22]])

    def test_coupled32_first_certified_power_and_exact_witnesses(self):
        result = direct_power_reference({
            "matrix": coupled(32), "tolerance": "1/100000000",
        })
        self.assertEqual(result["first_certified_power"], 11)
        self.assertEqual(
            result["bound_squared"], "505768790783/24414062500000000000000000000",
        )
        self.assertEqual(
            result["previous_bound_squared"], "217552976919/244140625000000000000000000",
        )
        self.assertTrue(all(
            Fraction(value) > Fraction(1, 10**16)
            for value in result["checked_bound_squared"][:-1]
        ))
        self.assertLessEqual(Fraction(result["bound_squared"]), Fraction(1, 10**16))

    def test_nonnormal8_control_matches_other_dimension(self):
        result = direct_power_reference({
            "matrix": coupled(8), "tolerance": "1/100000000",
        })
        self.assertEqual(result["first_certified_power"], 11)
        self.assertEqual(result["bound_squared"], "595333999/39062500000000000000000000")

    def test_nilpotent_geometric_failure_has_exact_power_success(self):
        result = direct_power_reference({"matrix": [[1,2],[0,1]]})
        self.assertEqual(result["first_certified_power"], 2)
        self.assertEqual(result["bound_squared"], "0")
        self.assertEqual(result["checked_bound_squared"], ["4","0"])

    def test_search_limit_and_invalid_maximum(self):
        result = direct_power_reference({
            "matrix": [["1/100"]], "tolerance": "1/100000000",
        }, max_power=64)
        self.assertEqual(result["status"], "reference_limit")
        self.assertIsNone(result["first_certified_power"])
        self.assertEqual(len(result["checked_bound_squared"]), 64)
        with self.assertRaises(CertificateError):
            direct_power_reference({"matrix": [[1]]}, max_power=True)

if __name__ == "__main__":
    unittest.main()
