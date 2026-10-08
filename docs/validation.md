# Validation

Version 0.4 was compiled with GNU Fortran using:

```text
-std=f2008 -Wall -Wextra -Werror -O2
```

All regression tests passed:

```text
homoclinic test: PASS
independent-networks test: PASS
nonuniform-grid test: PASS
single-saddle test: PASS
two-saddle heteroclinic test: PASS
```

## Original `tp2_0545` middle layer

The current rounded source field (`F10.4` precision) was processed independently by the v0.4 executable.
The program detected:

```text
saddles: 1
networks: 1
branches: 3
homoclinic: 1
heteroclinic: 0
boundary-going: 2
```

Refined saddle:

```text
x_s   =  9.3457212851E-04
y_s   = -2.1858213063E+00
psi_s = -1.0561806645E-01
det H = -4.4287131519E-04
```

Level-set quality:

```text
network RMS |psi-psi_s| = 1.632318E-06
network MAX |psi-psi_s| = 1.806645E-05
```

The generated point set was also compared with the previously generated reference separatrix point set:

```text
reference -> candidate RMS = 6.45142077E-03
candidate -> reference RMS = 4.43297882E-03
symmetric Hausdorff        = 6.29782239E-02
```

The Hausdorff maximum is below one x-grid spacing (0.08) of the rounded input field. The comparison should
be repeated after regenerating the source field with higher precision; near a hyperbolic point the topology
is especially sensitive to four-decimal rounding.

See `docs/images/tp2_0545_middle_validation.png` for the overlay.
