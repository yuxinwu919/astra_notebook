# AGENTS.md — Guidance for AI agents working on astra-notebook

astra-notebook is a Jupyter-first workbench for the DESY ASTRA particle
tracking code and its Generator. It replaces the official postpro /
lineplot / fieldplot graphics programs (not available on macOS) with a
modern frontend/backend split.

## Purpose and architecture

* **Frontend**: 6 task-specific Jupyter notebooks in `notebooks/`
  (01_generator, 02_astra, 03_postpro, 04_lineplot, 05_fieldplot,
  06_examples) - one per original program plus the examples summary.
  In addition, examples/<name>/<name>.ipynb holds one detailed teaching
  notebook per official example (8 total; shared spec in
  examples/_examples_spec.py).
  UI text and parameter help are in Chinese; plot labels are in English.
* **Backend**: the plain Python package `astra_tools/` (no packaging, no
  PyPI). The whole project folder is copied as-is to use it. Notebooks
  bootstrap the backend via `notebooks/_bootstrap.py` (locates the
  project root from `__file__`, so it works at any path).
* Some I/O code was vendored and adapted from lume-astra
  (ChristopherMayes/lume-astra, Apache-2.0); see module docstrings for
  attribution. Do NOT import lume-astra or pmd-beamphysics as runtime
  dependencies — the project is self-contained.

## Physics rules (audited against the ASTRA Manual V3.2 and validated
against real ASTRA output; DO NOT change without re-validating)

1. **Active particles**: `status > 1` (manual section 4.13). status 0/1
   = passive probes (excluded from statistics), -6..-1 = at cathode
   (not started), < -6 = lost.
2. **Distribution files**: the first row (ASCII) is the reference
   particle in ABSOLUTE coordinates; the remaining particles' z, pz and
   clock are RELATIVE to it and must be converted to absolute on read
   (see astra_tools/io/astra_dist.py). Binary files carry a 5-float
   header instead; pz/clock are still relative to it.
3. **Reference momentum**: header[1] / ref pz is MOMENTUM in eV/c, not
   kinetic energy. gamma = sqrt(1 + (p/mc)^2). Never use p/mc^2.
4. **Emittance in solenoid fields** (manual 4.13.1): use the canonical
   momentum p~x = px + c*Bz*y/2, p~y = py - c*Bz*x/2 with Bz = on-axis
   solenoid field at the bunch center. Divergences are formed with the
   reference momentum: x' = p~x / p_ref.
5. **Emittance units**: ASTRA prints "pi mm mrad"; numerically the
   value equals eps_n in mm.mrad (the pi marks the RMS phase-space
   ellipse AREA semantics: eps_rms = a*b, area = pi*a*b). Converting a
   file value to SI is a plain x1e-6 -> m.rad. NEVER multiply by pi.
   Display units in plots: "[pi mm mrad]" with values = eps_n*1e6.
   Longitudinal emittance: keV.mm (= eV.m), no pi factor numerically.
6. **Energy spread**: sigma_E/E from per-particle kinetic energies
   E_kin = sqrt(pz^2 + m^2 c^4) - m c^2, NOT sigma_p/p.
7. **Xemit corr column** = cov(u,u')/sigma_u (mrad), not the raw
   covariance. The Zemit counterpart <z E'>_avr is stored in keV and
   equals cov(z,E_kin)/sigma_z (x1e3 -> eV). Weighted moments use |q|
   weights so a mixed-sign bunch never silently falls back to
   unweighted statistics. The charge sign itself is kept internally
   and only dropped at the display/export boundary (|Q| shown).
8. All analysis returns SI internally; display-unit conversion happens
   only at the plot/export boundary.

## Testing discipline (five layers — all required for new code)

1. Unit tests for pure functions.
2. Golden-sample regression: archived outputs of the 9 official
   examples (see examples/); tests compare numbers, not requiring the
   ASTRA binaries.
3. Cross-validation vs ASTRA's own output: our statistics must match
   Xemit/Zemit columns within 0.5% (see test/test_cross_validation.py).
4. Plot correctness: axis labels/units/legends asserted
   (test/test_plots.py). Any new plot must assert its units.
5. End-to-end: every notebook must pass
   `jupyter nbconvert --to notebook --execute <nb> --ExecutePreprocessor.kernel_name=astra-notebook`
   (requires the ASTRA/Generator binaries).

Run all tests:  .venv/bin/python -m pytest test/ -q

## Conventions

* Python 3.9+ compatible syntax; numpy/scipy/matplotlib/pandas only for
  the core (ipywidgets for the frontend layer).
* Chinese docstrings/comments for user-facing modules are acceptable;
  plot labels and identifiers in English.
* Do not reintroduce packaging (pyproject build systems, console
  scripts, version pinning machinery). Copy-the-folder is the delivery
  model.
* Plot style: 2D phase-space plots are plain SCATTER plots with
  deterministic subsampling (max_points) and 0.5-99.5 percentile range
  clipping so outliers can never collapse the display; complete legends
  and unit strings everywhere. The RMS-ellipse overlay was removed
  (user decision 2026-08); the KDE density engine stays available in
  astra_tools/plot/_density.py for optional use, but is NOT the default.
* Git: inputs and small golden samples are committed; large regenerable
  outputs stay local (see .gitignore).

## Documentation

* docs/user_guide/  — Chinese user manual
* docs/dev_manual/   — Chinese developer manual
* docs/physics_notes/ — physics audit memos (formula -> manual citation
  -> numerical validation), one entry per audited convention
