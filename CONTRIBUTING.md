# Contributing

Contributions are welcome.

## Development workflow

1. Create a feature branch.
2. Keep the numerical core independent of any particular physical model.
3. Add or update a regression test for numerical/topological changes.
4. Run:

```bash
fpm test
```

5. For Python changes, run at least:

```bash
python -m py_compile python/*.py
```

## Numerical changes

Changes to saddle detection, contour reconstruction, graph assembly, or tolerance defaults
should include a short explanation of the numerical motivation in the pull request.

The test suite intentionally includes:

- a single saddle;
- a heteroclinic two-saddle network;
- two homoclinic loops;
- multiple independent separatrix networks.
