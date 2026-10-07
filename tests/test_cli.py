"""Exercise the real CLI with in-memory streams, without filesystem writes."""
import contextlib
import io
import json
import unittest
from unittest.mock import patch
from nrcert import make_certificate
from nrcert.checker import check_json
from nrcert.cli import main

def invoke(arguments, input_text):
    output, errors = io.StringIO(), io.StringIO()
    with patch("sys.stdin", io.StringIO(input_text)):
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            code = main(arguments)
    return code, output.getvalue(), errors.getvalue()

class CommandLineCertificates(unittest.TestCase):
    def test_certify_then_check_stdin(self):
        code, output, errors = invoke(["certify", "-"], '{"matrix":[[1]]}')
        self.assertEqual(code, 0)
        self.assertEqual(errors, "")
        self.assertEqual(check_json(output)["target_certified"], True)
        code, checked, errors = invoke(["check", "-"], output)
        self.assertEqual(code, 0)
        self.assertTrue(json.loads(checked)["valid"])

    def test_conditional_flag_required_on_both_commands(self):
        request = '{"kind":"polynomial","matrix":[[1]],"coefficients":[1],"mode":"conditional_two"}'
        self.assertEqual(invoke(["certify", "-"], request)[0], 1)
        code, output, errors = invoke(["certify", "-", "--allow-conditional"], request)
        self.assertEqual(code, 0)
        self.assertEqual(invoke(["check", "-"], output)[0], 1)
        self.assertEqual(invoke(["check", "-", "--allow-conditional"], output)[0], 0)

    def test_failure_returns_two_and_contains_evidence(self):
        code, output, errors = invoke(["certify", "-"], '{"matrix":[[1,2],[0,1]]}')
        self.assertEqual(code, 2)
        self.assertEqual(json.loads(output)["evidence"]["status"], "no_geometric_contraction")

    def test_bad_json_and_wrong_matrix_fail(self):
        for text in ("{", '{"matrix":[[1.1]]}', "[]"):
            with self.subTest(text=text):
                code, output, errors = invoke(["certify", "-"], text)
                self.assertEqual(code, 1)
                self.assertEqual(output, "")
                self.assertIn("nrcert: error:", errors)

    def test_json_consumer_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            check_json("{")

if __name__ == "__main__":
    unittest.main()

