# FEM Learning Laboratory (monochrome web application)

A simple **single-server** educational FEM application. Enter nodes, members, supports, nodal loads, and stiffness values. Click **SOLVE STRUCTURE** to solve an indeterminate bar/truss/beam/frame using the Python backend. Click **SHOW COMPLETE SOLUTION** to inspect the real shape functions, element stiffness Kₑ, transformations, assembled global Kᴳ, boundary-condition matrix H, Lagrange multiplier augmented system, solved displacements and reactions.

**No Node.js, React, npm, or separate frontend server is required.**

Use Python 3.12 or later. Dependencies are declared once in `pyproject.toml`.

## Windows: quick start

**Easiest method:** Double-click `START_FEM_LAB.bat` inside the extracted project folder. The script creates a virtual environment (if necessary), installs dependencies, starts the server on **port 8010**, and opens your browser. Keep the black terminal window open. Port 8010 is separate from BeamFEM Studio’s port 8000.

**Manual method:** Open **one PowerShell terminal** in VS Code, from the extracted `FEM_Learning_Laboratory` folder. Run **only** the commands, never the `PS D:\...>` prompt label.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m uvicorn backend.app:app --reload --port 8010
```

Open **http://127.0.0.1:8010** in Chrome or Edge. API docs at **http://127.0.0.1:8010/docs**. Leave the terminal running.

For later starts, just run:

```powershell
.\.venv\Scripts\python.exe -m uvicorn backend.app:app --reload --port 8010
```

Linux / macOS:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -m uvicorn backend.app:app --reload --port 8010
```

## What is implemented

| Structural model | Physical node DOFs | User-input minimum | p-refinement | h-refinement |
|---|---|---|---|---|
| 1D axial bar | ux | EA | 1–6, genuine C0 axial bubble modes | 1–8 subdivisions per member |
| 2D truss | ux, uy | EA | 1–6, genuine C0 axial bubble modes | 1–8 |
| 3D truss | ux, uy, uz | EA | 1–6, genuine C0 axial bubble modes | 1–8 |
| 1D Euler–Bernoulli beam | uy, rz | EI | 3–6, Hermite with genuine C1 bending bubbles | 1–8 |
| 2D Euler–Bernoulli frame | ux, uy, rz | EA, EI | 3–6, axial C0 + bending C1 | 1–8 |
| 3D Euler–Bernoulli space frame | ux, uy, uz, rx, ry, rz | EA, EIy, EIz, GJ | 3–6, axial/torsion C0 + two bending C1 families | 1–8 |

- Nodes and element inputs are fully editable, including individual p and h refinement and fixed DOFs.
- Supports: free, fully fixed, pinned/roller presets, or individually restrained DOFs.
- Nodal loads/forces and moments; local uniformly distributed member loads for beam/frame formulations.
- Static linear analysis using the **Lagrange multiplier method** `[Kᴳ Hᵀ; H 0][d; λ]=[F; g]`, solved with NumPy/SciPy and numerical equilibration; SymPy/native MathML step-by-step solution viewer.
- Local element stiffness, transformations, local/global member forces, full global stiffness assembly, explicit H matrix, simultaneous displacement/multiplier solve, and support reactions from `R=Kᴳd−F=−Hᵀλ`.
- Built-in stable example for each of six model types.
- Save/load models as JSON; download a reproducible Python script; print the manual.
- Black-and-white academic visual design. Chrome 156+ and modern Edge render the equations using **native MathML** (no internet/CDN required).

### How to use

1. Choose **Structure type** and press **Load example**.
2. Change nodes x/y/z, support types and nodal forces/moments.
3. Change each member's node connectivity i/j, EA/EI (or EIy/EIz/GJ), p and h; use qy for beam distributed loads.
4. Click **SOLVE STRUCTURE** to display displacements, reactions, end forces and a schematic.
5. Click **SHOW COMPLETE SOLUTION**; inspect 10 steps in this order: shape functions → element stiffness Kₑ → transformation T → global stiffness Kᴳ → constraint matrix H → Lagrangian → augmented matrix → displacements and multipliers → reactions → member forces.
6. If you choose an excessively flexible or unsupported structure, the solver will explain why its system is singular.

### Theory / limitations

- Base units are SI: m, N, rad; EA in N, EI and GJ in N·m², distributed load in N/m.
- Material/section values enter as stiffness products directly, e.g. **EA = E·A**; 1D/2D/3D frame bending requires **EI**, 3D frame requires **EIy, EIz and GJ**. EA alone cannot describe frame bending/torsion.
- 1D beam orientation is along global X. 2D members lie in XY. 3D coordinates may have arbitrary orientation. The 3D drawing is an oblique projection; edit numeric coordinates for accuracy.
- Distributed member loads are in **local axes**; axial/truss examples accept concentrated nodal loads only. This is not a nonlinear, dynamic, buckling or design-code checker.
- h-refinement physically subdivides each element; p-refinement changes its compatible interpolation basis and introduces additional internal DOFs. The internal modes are solved explicitly. **Nodal-load-only bar/truss solutions can be identical under p refinement when the exact displacement field is already linear**; that is expected, not an error.
- Stiffness matrices above 26 DOFs are suppressed in the manual for readability, while still solved numerically. Solver demonstration limit: 550 DOFs. Visual deformations are exaggerated, not to engineering scale.
- This prototype is for education and must be independently validated before being used for safety-critical engineering design. It has numerical benchmark tests but no third-party software certification.

## Run tests

```powershell
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m pytest -q
```

On Linux / macOS:

```bash
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m pytest -q
```

`pip install -e .` installs this project and its runtime dependencies. The `-e`
option keeps it linked to the source folder, where the app also reads `static/`.
`.[test]` additionally installs the optional test tools. Add runtime packages to
`[project].dependencies`; add test tools to `[project.optional-dependencies].test`.
`requirements.txt` is a compatibility pointer to this project, not a second list.

## Deploy to Vercel

Use the FastAPI framework preset with the repository root as the Root Directory.
`[tool.vercel]` in `pyproject.toml` points to `backend.app:app`, while `[project]`
declares the runtime dependencies. No custom build command is needed.
After deployment, check `/api/health`, the homepage, and solving a built-in example.

## Main files

- `backend/models.py` — input model/schema validation.
- `backend/elements.py` — compatible polynomial/Hermite interpolation and 1D/2D/3D elements.
- `backend/solver.py` — h mesh expansion, global assembly, explicit H constraints, scaled augmented Lagrange multiplier solver, reactions.
- `backend/derivations.py` — math solution manual using presentation MathML.
- `backend/app.py` — FastAPI endpoints and website static serving.
- `static/index.html`, `static/style.css`, `static/app.js` — dependency-free black-and-white web interface.
- `tests/test_solver.py` — analytical and regression tests.

## Version 0.2 — Lagrange multiplier correction

The original version used `K_ff d_f=F_f` (partitioning) and described it in the solution manual. This version **actually solves** the augmented system instead:

```text
| K_G  H^T | | d      |   | F |
| H    0  | | lambda | = | g |
```

- `H` has one independent row per restrained displacement DOF, with a `1` in that global DOF's column.
- `g=0` for the current fixed/roller/pin support inputs; prescribed nonzero settlements and general multipoint constraints are not yet implemented.
- The energy convention is `Pi(d,lambda) = (1/2)d^T K_G d - d^T F + lambda^T(Hd-g)`, so the physical support reaction is **`R=K_G d-F=-H^T lambda`**. Multipliers are therefore the **negative** of the reaction for a unit-row support constraint.
- The implementation uses reversible diagonal equilibration internally, but the solution manual shows the **unscaled physical** K_G, H and augmented matrices.
- The 10-step manual shows the actual selected element matrices, assembled global matrices, constraints, displacements, multipliers, and reaction checks.
- Reversed 1D beam connectivity now transforms nodal rotation signs correctly.
- Automated numerical/API tests include analytical cantilever UDL, indeterminate bars, mixed p/h refinement, matrix equilibrium and multiplier reaction signs.

**Updating an existing installation:** Stop the old server with Ctrl+C, extract the *update ZIP* directly into your existing `FEM_Learning_Laboratory` folder, allow Windows to replace the matching files, then double-click `START_FEM_LAB.bat`. Existing `.venv` and saved model JSON files can remain.
