# Nonautonomous invariant-manifold extension

This document defines the first prototype architecture for extending Separatrix Toolkit from stationary level-set separatrices to stable/unstable manifolds of a hyperbolic periodic point of a 2-D Poincare map.

The stationary solver remains unchanged. The new layer works with a generic map callback

\[
P:\;\mathbf{x}_0\mapsto\mathbf{x}_1,
\]

so a hydrodynamic model can provide the map by integrating one forcing period with its own ODE solver.

Prototype sequence:

1. refine a periodic point from `P(x)-x=0`;
2. estimate the monodromy matrix `DP` by centered finite differences;
3. classify the point as hyperbolic from the eigenvalues of `DP`;
4. recover stable and unstable eigendirections;
5. seed short material segments along those directions;
6. iterate unstable seeds forward and stable seeds backward;
7. detect polyline intersections;
8. form a candidate turnstile lobe and compute its signed polygon area.

The first implementation intentionally keeps the Poincare map abstract. This avoids coupling the reusable toolkit to a particular oceanographic model or to DOP853. A model-specific adapter can be supplied by a caller.

## Numerical cautions

- A time-periodic flow requires a map inverse for backward stable-manifold continuation. If no analytic inverse is available, the adapter should integrate the ODE backward by one forcing period.
- Finite-difference Jacobians should be checked against at least two perturbation sizes.
- Manifold iteration needs adaptive point insertion once adjacent points become too far apart. The initial prototype does not yet claim production-quality continuation.
- Lobe selection is a topological task: not every consecutive pair of manifold intersections defines the physically relevant turnstile lobe.
- The lobe-area predictor used in transport modelling must be computed independently of the transport ensemble to avoid outcome leakage.

The first target application is the regime `layer=1`, `omega0=0.024`, `epsilon=0.05`.
