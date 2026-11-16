"""
main.py - The browser side of the website (runs inside PyScript).

1. read the numbers typed into the page,
2. call the fem package (the same code that pytest checks),
3. put the results and figures back on the page.
"""

import asyncio
import base64
import io

import matplotlib.pyplot as plt
from pyscript import document, when

from fem import PRESETS, load_sweep, run_case
from fem.plotting import plot_field, plot_load_sweep, plot_top_stress

MAX_DOF = 2000            # keeps the dense solve fast in the browser
last = {}                 # most recent result, so the field menu can redraw without re-solving


def el(id_):
    return document.getElementById(id_)


def num(id_):
    try:
        return float(el(id_).value)
    except ValueError:
        raise ValueError(f"'{id_}' must be a number")


def set_status(text, kind="info"):
    el("status").innerText = text
    el("status").className = f"status {kind}"


def show_figure(fig, img_id):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    img = el(img_id)
    img.src = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()
    img.style.display = "block"


def read_inputs():
    nx, ny = int(num("nx")), int(num("ny"))
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be at least 1")
    if 2 * (nx + 1) * (ny + 1) > MAX_DOF:
        raise ValueError(f"Mesh too large for the browser (limit {MAX_DOF} DOF). Reduce nx or ny.")
    return dict(L=num("L"), h=num("h"), t=num("t"), p=num("p") * 1e3, load=el("load").value,
                E=num("E") * 1e9, nu=num("nu"), yield_strength=num("sy") * 1e6,
                nx=nx, ny=ny, element=el("element").value)


def results_table(r):
    bending = r["inputs"]["load"] == "bending"
    name = "Tip deflection" if bending else "Elongation"
    rows = [
        ("Mesh", f"{r['n_elements']} {r['inputs']['element']} elements, {r['n_dof']} DOF"),
        (f"{name}, FEM", f"{r['tip']*1e3:.4g} mm"),
        (f"{name}, beam theory", f"{r['tip_ref']*1e3:.4g} mm ({r['tip_error_pct']:+.2f} %)"),
        ("Max von Mises stress", f"{r['max_von_mises']/1e6:.4g} MPa"),
        ("Max strain |εxx|", f"{r['max_strain']*1e6:.4g} µε"),
        ("Safety factor vs yield", f"{r['safety_factor']:.2f}"
         + ("  ⚠ yields" if r["safety_factor"] < 1 else "")),
    ]
    if bending:
        rows.insert(4, ("σxx top, mid-span: FEM / theory",
                        f"{r['sigma_mid']/1e6:.4g} / {r['sigma_mid_ref']/1e6:.4g} MPa"))
    rows.append(("Equilibrium check: Σ reactions",
                 f"Rx = {r['reaction_x']:.4g} N, Ry = {r['reaction_y']:.4g} N"))
    el("results").innerHTML = "<table>" + "".join(
        f"<tr><th>{a}</th><td>{b}</td></tr>" for a, b in rows) + "</table>"


def apply_preset():
    key = el("preset").value
    custom = key == "custom"
    for f in ("E", "nu", "sy"):
        el(f).readOnly = not custom
    if not custom:
        m = PRESETS[key]
        el("E").value = f"{m['E']/1e9:.4g}"
        el("nu").value = f"{m['nu']:.3g}"
        el("sy").value = f"{m['yield']/1e6:.4g}"


def draw_field():
    if last:
        show_figure(plot_field(last["r"], el("field").value), "img-field")


async def solve():
    set_status("Solving…", "busy")
    el("solve-btn").disabled = True
    await asyncio.sleep(0.05)                 # let the browser repaint first
    try:
        r = run_case(**read_inputs())
        last["r"] = r
        results_table(r)
        draw_field()
        show_figure(plot_load_sweep(load_sweep(r), r["yield"]), "img-sweep")
        if r["inputs"]["load"] == "bending":
            show_figure(plot_top_stress(r), "img-top")
        else:
            el("img-top").style.display = "none"
        set_status(f"Solved {r['n_dof']} DOF in {r['solve_time_s']*1e3:.0f} ms.", "ok")
    except Exception as exc:                  # show any problem on the page
        set_status(f"Error: {exc}", "error")
    el("solve-btn").disabled = False


# ----------------------------------------------------------------- events
# PyScript's @when passes the browser event to the handler, so handlers take
# `event` and call the plain functions above (which start-up code can also call).
@when("change", "#preset")
def on_preset(event):
    apply_preset()


@when("change", "#field")
def on_field(event):
    draw_field()


@when("click", "#solve-btn")
async def on_solve(event):
    await solve()


# ----------------------------------------------------------------- start-up
apply_preset()
el("loading").style.display = "none"
el("solve-btn").disabled = False
asyncio.ensure_future(solve())
