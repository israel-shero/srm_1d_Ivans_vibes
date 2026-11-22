"""Mocked tests for claim notice publication; no GitHub writes."""
import importlib.util
import pathlib
import sys
import unittest
from unittest.mock import patch

HERE = pathlib.Path(__file__).parent
sys.path.insert(0,str(HERE))
import publish_claim as mod

class PublishTests(unittest.TestCase):
    def args(self):
        return ["owner/repo","INF-01","codex-A","work/INF-01/codex-A","Profile functions","human-operator","1234567","a"*40,6,"token"]
    def test_notice_contains_workflow_and_owner(self):
        with patch.object(mod,"request",return_value=(201,{"id":123,"html_url":"https://github.com/example"})) as req:
            r=mod.publish(*self.args())
            self.assertEqual(r["status"],"PUBLISHED")
            self.assertEqual(req.call_args.args[0],"POST")
            self.assertIn("issues/6/comments",req.call_args.args[1])
            b=req.call_args.args[3]["body"]
            self.assertIn("codex-A",b);self.assertIn("actions/runs/1234567",b)
            self.assertIn("human-operator",b)
    def test_incomplete_response_fails(self):
        with patch.object(mod,"request",return_value=(201,{})):
            with self.assertRaises(mod.ClaimError):mod.publish(*self.args())
    def test_http_failure_fails(self):
        with patch.object(mod,"request",side_effect=mod.ClaimError("403")):
            with self.assertRaises(mod.ClaimError):mod.publish(*self.args())
    def test_bad_branch_fails_without_network(self):
        a=self.args();a[3]="../unsafe"
        with patch.object(mod,"request") as req:
            with self.assertRaises(mod.ClaimError):mod.publish(*a)
            req.assert_not_called()
    def test_scope_newline_rejected(self):
        a=self.args();a[4]="hello\nmalicious"
        with self.assertRaises(mod.ClaimError):mod.publish(*a)
    def test_agent_label_rejected(self):
        a=self.args();a[2]="agent with spaces"
        with self.assertRaises(mod.ClaimError):mod.publish(*a)
if __name__=="__main__":
    unittest.main()
