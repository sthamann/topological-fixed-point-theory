import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from semantic_explorer import review_semantics as review


def packet():
    return {
        "claim_id": "claim:EM.FP.01",
        "semantic": {"interpretation_de": "sample"},
        "source_excerpts": [{"source_id": "source:a", "path": "source.py", "sha256": "digest",
            "sha256_verified": True, "text": "A short exact source quote."}],
        "unavailable_sources": [],
    }


def clean_response():
    return {
        "schema_version": 1,
        "claim_id": "claim:EM.FP.01",
        "review_state": "no_findings",
        "reviewed": False,
        "issues": [],
        "coverage": {
            "fields_checked": ["interpretation_de"],
            "source_ids_checked": ["source:a"],
            "unavailable_sources": [],
            "limitations": [],
        },
        "notice": "Unreviewed second-pass proposal; not an authoritative correction or proof verdict.",
    }


def write_graph(path, node, sources=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(__import__("json").dumps({"nodes": [node], "edges": [], "sources": sources or []}), encoding="utf-8")


class ReviewSemanticsTests(unittest.TestCase):
    def test_successful_no_findings_review_is_proposal_only(self):
        with patch.object(review.extract, "resolve_key", return_value="secret"), patch.object(
            review.extract, "api", return_value={"status": "completed", "output": [{"content": [{"type": "output_text", "text": __import__("json").dumps(clean_response())}]}]}
        ) as api:
            result = review.run_review(packet(), "test-model", "/unused")
        self.assertTrue(result["proposal_only"])
        self.assertEqual(result["review"]["review_state"], "no_findings")
        self.assertEqual(api.call_args.args[0], "responses")

    def test_invalid_quote_fails_closed_without_output(self):
        response = clean_response()
        response["review_state"] = "reviewed_with_findings"
        response["issues"] = [{
            "field_path": "interpretation_de", "code": "source_contradiction", "severity": "high",
            "finding": "Mismatch", "suggested_correction": "Recheck", "evidence": [
                {"source_id": "source:a", "relation": "contradicts", "quote": "not in source"}
            ]
        }]
        with patch.object(review.extract, "resolve_key", return_value="secret"), patch.object(
            review.extract, "api", return_value={"status": "completed", "output": [{"content": [{"type": "output_text", "text": __import__("json").dumps(response)}]}]}
        ):
            with tempfile.TemporaryDirectory() as directory:
                out = Path(directory) / "proposal.json"
                with self.assertRaisesRegex(ValueError, "validation failed"):
                    proposal = review.run_review(packet(), "test-model")
                    review.atomic_write_new(out, proposal)
                self.assertFalse(out.exists())

    def test_no_findings_rejected_when_source_is_unavailable(self):
        bad_packet = packet()
        bad_packet["unavailable_sources"] = [{"source_id": "source:missing", "reason": "unavailable"}]
        self.assertTrue(any("no_findings is not allowed" in error for error in review.validate_response(clean_response(), bad_packet)))

    def test_empty_annotation_preflight_prevents_key_resolution_and_api_call(self):
        bad_packet = {"claim_id": "claim:FLAV.CKM.01", "semantic": {}, "source_excerpts": []}
        with patch.object(review.extract, "resolve_key") as resolve, patch.object(review.extract, "api") as api:
            with self.assertRaisesRegex(ValueError, "no populated semantic annotation"):
                review.run_review(bad_packet, "test-model", "/unused")
        resolve.assert_not_called()
        api.assert_not_called()

    def test_unavailable_sources_fail_preflight_before_key_resolution(self):
        bad_packet = packet()
        bad_packet["unavailable_sources"] = [{"source_id": "source:missing", "reason": "missing"}]
        with patch.object(review.extract, "resolve_key") as resolve, patch.object(review.extract, "api") as api:
            with self.assertRaisesRegex(ValueError, "referenced sources are unavailable"):
                review.run_review(bad_packet, "test-model", "/unused")
        resolve.assert_not_called()
        api.assert_not_called()

    def test_priority_pending_claim_falls_back_to_enriched_base_and_location_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root / "semantic_explorer" / "data"
            priority_path = data / "priority" / "knowledge.json"
            base_path = data / "knowledge.json"
            repo = root / "repo"
            (repo / "verification").mkdir(parents=True)
            (repo / "verification" / "status_ledger.csv").write_text("claim,CKM core\n", encoding="utf-8")
            (repo / "tfpt_1_architecture.tex").write_text(
                "intro\nCKM parameters\n\\delta_{\\mathrm{CKM}}=\\pi/3+3\\lC^2\n"
                "$s_{12}^{\\mathrm{CKM}}$ = \\lambda_C\ncontext\n", encoding="utf-8")
            ledger = {"claim_id": "FLAV.CKM.01", "claim": "CKM core delta=pi/3+3lC^2",
                "location": "tfpt_1", "script": "-"}
            source = {"id": "source:ledger", "path": "verification/status_ledger.csv",
                "sha256": review.sha256((repo / "verification" / "status_ledger.csv").read_bytes()),
                "line_start": 1, "line_end": 1}
            semantic = {"title_de": "CKM phase", "overview_de": "Numerischer Readout",
                "interpretation_de": "Keine unabhängige Herleitung.", "results": [],
                "proof_gaps": [], "source_refs": ["source:ledger"]}
            write_graph(priority_path, {"id": "claim:FLAV.CKM.01", "type": "claim", "attrs": {
                "enrichment_state": "pending", "semantic": {}, "ledger": ledger}})
            write_graph(base_path, {"id": "claim:FLAV.CKM.01", "type": "claim", "attrs": {
                "enrichment_state": "enriched", "semantic": semantic, "ledger": ledger,
                "provenance": {"created_at": "2026-10-12T00:00:00Z"}}}, [source])
            packet = review.build_packet(priority_path, repo, "FLAV.CKM.01")
            self.assertEqual(packet["semantic"]["title_de"], "CKM phase")
            self.assertEqual(packet["source_annotation_sha256"], review.sha256(
                __import__("json").dumps(semantic, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()))
            self.assertTrue(any(x["path"] == "tfpt_1_architecture.tex" for x in packet["source_excerpts"]))
            self.assertTrue(all(x["sha256_verified"] for x in packet["source_excerpts"]))
            excerpts = "\n".join(x["text"] for x in packet["source_excerpts"])
            self.assertIn("delta_{\\mathrm{CKM}}", excerpts)
            self.assertIn("s_{12}^{\\mathrm{CKM}}", excerpts)
            tex = next(x for x in packet["source_excerpts"] if x["path"] == "tfpt_1_architecture.tex")
            tex_hash = review.sha256((repo / "tfpt_1_architecture.tex").read_bytes())
            expected_source_id = "source:" + review.sha256(
                f"tfpt_1_architecture.tex:{tex_hash}:{tex['line_start']}:{tex['line_end']}".encode())[:20]
            self.assertEqual(tex["source_id"], expected_source_id)

    def test_malformed_payload_fails_closed_without_validator_crash(self):
        response = clean_response()
        response["coverage"]["fields_checked"] = [["not", "a", "string"]]
        response["issues"] = {"unexpected": "object"}
        errors = review.validate_response(response, packet())
        self.assertTrue(any("arrays of strings" in error for error in errors))
        self.assertTrue(any("issues must be an array" in error for error in errors))

    def test_atomic_output_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "proposal.json"
            out.write_text("existing", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                review.atomic_write_new(out, {"proposal_only": True})
            self.assertEqual(out.read_text(encoding="utf-8"), "existing")


if __name__ == "__main__":
    unittest.main()
