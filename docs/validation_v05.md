# Validation

Version 0.5 was checked with `gfortran` using:

```text
-std=f2008 -Wall -Wextra -Werror -O2
```

All Fortran regression tests passed:

```text
homoclinic test: PASS
independent-networks test: PASS
nonuniform-grid test: PASS
single-saddle test: PASS
two-saddle heteroclinic test: PASS
```

The installable Python package was also smoke-tested locally with:

```bash
pip install -e . --no-build-isolation
sepviz inspect <separatrices.dat>
```

For the hydrodynamic middle-layer validation dataset the Python inspector reported:

```text
layers: 1
networks: 1
branches: 3
points: 490
```

The corresponding publication-style figure is stored in:

```text
docs/images/v05_middle_demo.png
```

The numerical core remains independent of this hydrodynamic example.
