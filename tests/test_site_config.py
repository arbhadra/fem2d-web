"""Checks that the browser configuration loads every solver module (with the same cache-busting version)."""

import glob
import os
import re

ROOT = os.path.join(os.path.dirname(__file__), "..")


def test_pyscript_toml_lists_every_fem_module_with_one_version():
    toml = open(os.path.join(ROOT, "pyscript.toml")).read()
    html = open(os.path.join(ROOT, "index.html")).read()
    modules = {os.path.basename(p) for p in glob.glob(os.path.join(ROOT, "fem", "*.py"))}
    listed = set(re.findall(r'"\./fem/([a-z_]+\.py)\?v=', toml))
    assert modules == listed, f"missing from pyscript.toml: {modules - listed}"
    versions = set(re.findall(r"\?v=([0-9.]+)", toml + html))
    assert len(versions) == 1, f"mixed versions: {versions}"
