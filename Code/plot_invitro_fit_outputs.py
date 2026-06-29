#!/usr/bin/env python3
"""Compatibility entrypoint for the in vitro PKPD/live-dead plotter."""

from pathlib import Path
import runpy


MODULE_SCRIPT = (
    Path(__file__).resolve().parent
    / "in-vitro"
    / "pkpd_live_dead_model"
    / "plot_invitro_fit_outputs.py"
)

runpy.run_path(str(MODULE_SCRIPT), run_name="__main__")
