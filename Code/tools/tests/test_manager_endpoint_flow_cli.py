from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
TOOLS_DIR = REPO_ROOT / "Code/tools"
sys.path.insert(0, str(TOOLS_DIR))

from materialize_figure_assets import PANEL_SPECS  # noqa: E402


class ManagerEndpointFlowCliTest(unittest.TestCase):
    def run_manager(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(REPO_ROOT / "Manager.sh"), *args],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
        )

    def test_endpoint_flow_is_opt_in(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text(encoding="utf-8")
        default_line = next(
            line for line in manager_text.splitlines() if line.startswith('modules="')
        )
        self.assertNotIn("in_vivo_endpoint_flow", default_line)

    def test_endpoint_flow_materializes_only_the_two_panel_composite(self) -> None:
        specs = [
            spec
            for spec in PANEL_SPECS
            if spec["module"] == "in_vivo_endpoint_flow"
        ]
        self.assertEqual(
            [spec["panel"] for spec in specs],
            ["SuppFig9", "SuppFig9_png"],
        )
        self.assertEqual(
            {spec["source"] for spec in specs},
            {
                "figures/panel_SuppFig9_endpoint_flow_cytometry.pdf",
                "figures/panel_SuppFig9_endpoint_flow_cytometry.png",
            },
        )
        manager_text = (REPO_ROOT / "Manager.sh").read_text(encoding="utf-8")
        panels_only_case = manager_text.split(
            'if [[ "${mode}" == "panels-only" ]]', 1
        )[1].split('if [[ "${skip_analysis_loop}"', 1)[0]
        self.assertIn("in_vivo_endpoint_flow", panels_only_case)

    def test_check_only_resolves_portable_reviewed_inputs(self) -> None:
        result = self.run_manager(
            "--mode",
            "check-only",
            "--modules",
            "in_vivo_endpoint_flow",
            "--run-id",
            "endpoint_flow_contract",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(
            "Code/in-vivo/flow_cytometry/run_endpoint_flow.sh", result.stdout
        )
        wrapper_text = (
            REPO_ROOT / "Code/in-vivo/flow_cytometry/run_endpoint_flow.sh"
        ).read_text(encoding="utf-8")
        self.assertIn("extract_endpoint_flow.py", wrapper_text)
        self.assertIn("reconstruct_endpoint_flow.R", wrapper_text)
        manager_text = (REPO_ROOT / "Manager.sh").read_text(encoding="utf-8")
        self.assertIn("scripts/agentRrunner.sh", manager_text)
        self.assertIn(
            "Data/in-vivo/flow_cytometry/endpoint_tumors_20250128/crosswalk.tsv",
            result.stdout,
        )
        self.assertIn(
            "Data/in-vivo/flow_cytometry/endpoint_tumors_20250128/content_manifest.tsv",
            result.stdout,
        )
        self.assertIn(
            "Data/in-vivo/flow_cytometry/endpoint_tumors_20250128/paired_acquisition_sensitivity.tsv",
            result.stdout,
        )
        self.assertIn("--expected-samples 16", result.stdout)
        self.assertIn("--min-human-cells 1000", result.stdout)
        self.assertNotIn("/Volumes/Flow_Cytometry", result.stdout)


if __name__ == "__main__":
    unittest.main()
