"""
Runs main.py (the website code) against a fake browser page, so start-up and
button-handler errors are caught by CI before they reach the live site.
"""

import asyncio
import os
import runpy
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "fake_pyscript"))


def test_website_starts_solves_and_handles_events():
    import pyscript

    async def run():
        runpy.run_path(os.path.join(ROOT, "main.py"))
        await asyncio.sleep(0.5)                      # let the start-up solve finish
        page = pyscript.document.getElementById
        assert page("status").className == "status ok", page("status").innerText
        assert page("solve-btn").disabled is False
        assert page("img-field").src.startswith("data:image/png")

        page("field").value = "exx"
        pyscript.handlers["#field"]({})               # real PyScript passes an event
        page("preset").value = "custom"
        pyscript.handlers["#preset"]({})
        assert page("E").readOnly is False
        page("load").value = "tension"
        await pyscript.handlers["#solve-btn"]({})
        assert page("status").className == "status ok", page("status").innerText

    cwd = os.getcwd()
    os.chdir(ROOT)
    try:
        asyncio.run(run())
    finally:
        os.chdir(cwd)
