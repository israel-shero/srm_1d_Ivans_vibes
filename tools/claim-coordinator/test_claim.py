"""Offline tests. No network or GitHub permissions required."""
import importlib.util
import pathlib
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("claim", pathlib.Path(__file__).with_name("claim.py"))
claim = importlib.util.module_from_spec(spec)
spec.loader.exec_module(claim)

class ClaimTests(unittest.TestCase):
    def test_valid_claim(self):
        sha = "a"*40
        with patch.object(claim,"request",return_value=(201,{"ref":"refs/heads/claims/UQ-01","object":{"sha":sha}})) as call:
            r=claim.run("owner/repo","UQ-01",sha,"token","https://api.github.com","acquire")
            self.assertEqual(r["status"],"ACQUIRED")
            self.assertEqual(call.call_args.args[0],"POST")
            self.assertEqual(call.call_args.args[3],{"ref":"refs/heads/claims/UQ-01","sha":sha})
    def test_conflict_fails_closed(self):
        with patch.object(claim,"request",side_effect=claim.ClaimError("GitHub HTTP 422")):
            with self.assertRaises(claim.ClaimError):
                claim.run("owner/repo","UQ-01","a"*40,"token","https://api.github.com","acquire")
    def test_ambiguous_fails_closed(self):
        with patch.object(claim,"request",side_effect=claim.ClaimError("timeout")):
            with self.assertRaises(claim.ClaimError):
                claim.run("owner/repo","UQ-01","a"*40,"token","https://api.github.com","acquire")
    def test_invalid_task_rejected(self):
        with self.assertRaises(claim.ClaimError):
            claim.run("owner/repo","../escape","a"*40,"token","https://api.github.com","acquire")
    def test_no_base_rejected(self):
        with self.assertRaises(claim.ClaimError):
            claim.run("owner/repo","UQ-01",None,"token","https://api.github.com","acquire")
    def test_wrong_response_rejected(self):
        with patch.object(claim,"request",return_value=(201,{"ref":"refs/heads/claims/OTHER","object":{"sha":"a"*40}})):
            with self.assertRaises(claim.ClaimError):
                claim.run("owner/repo","UQ-01","a"*40,"token","https://api.github.com","acquire")
if __name__=="__main__":
    unittest.main()
