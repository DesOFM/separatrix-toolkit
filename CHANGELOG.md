# Changelog

## 0.8.0

- selected Apache License 2.0 for the public release;
- added `NOTICE` and final author/repository metadata;

- added a reproducible grid-convergence study;
- added shifted Duffing homoclinic and two-saddle heteroclinic analytic benchmarks;
- documented approximately second-order far-field geometric convergence;
- documented the local saddle neighborhood as the current uniform-geometry accuracy limit;
- added committed convergence CSV/JSON results and convergence plots;
- added `tools/run_convergence.py` and optional `convergence` Python dependency.


## 0.7.0

- added `REFERENCES.bib` with the 1991 and 1997 Sokolovskiy model references;
- added `docs/validation_model.md` explaining model provenance;
- clarified in README that the hydrodynamic model is cited as a validation-field source,
  while the separatrix algorithm is an independent implementation;
- kept scientific references separate from the software license.


## 0.6.0

- README restructure;
- clearer statement of scientific scope and limitations;
- GitHub repository positioning guide and suggested topics;
- issue templates and pull-request checklist;
- expanded Python package metadata;
- wording changed from "arbitrary" to "multiple" where appropriate to avoid overclaiming;
- explicit separation between the reusable numerical core and motivating research data.

## 0.5.0

- added installable Python package `separatrix-viz`;
- added unified `sepviz` CLI with `plot`, `prepare`, `compare`, and `inspect`;
- added PNG/PDF/SVG export;
- added optional background streamfunction contours and saddle annotations;
- added a quick-start Jupyter notebook;
- added `CITATION.cff` and gallery documentation;
- expanded CI with a complete Fortran-to-Python end-to-end smoke test.

## 0.4.0

- added spatial hashing for contour vertex and duplicate-edge assembly;
- added rectilinear nonuniform-grid support in numerical differentiation and scale selection;
- added per-branch arc length and level-set residual quality metrics;
- added network-level RMS/max residual diagnostics;
- added nonuniform-grid regression test;
- added Python reference-curve comparison utility with RMS and Hausdorff metrics.

## 0.3.0

- added command-line configuration;
- added optional Fortran namelist configuration;
- added automatic diagnostics report;
- added topology classification for homoclinic, heteroclinic and boundary-going branches;
- added homoclinic Duffing regression test;
- added multiple-independent-networks regression test;
- added GitHub Actions CI;
- added license file and contributing guide;
- documented the numerical algorithm and validation scope;
- retained a single `separatrices.dat` output format with no hard-coded topology.