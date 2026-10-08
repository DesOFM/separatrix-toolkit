# Numerical algorithm

The library operates on a scalar field \(\psi(x,y)\) sampled on a regular rectangular grid.

## 1. Critical-point detection

The discrete gradient is evaluated on the grid. Cells that bracket simultaneous zeros of

\[
\psi_x = 0,\qquad \psi_y = 0
\]

are treated as critical-point candidates. A bilinear approximation of the two gradient
components is solved inside each candidate cell.

## 2. Saddle refinement

Each candidate is refined with a local quadratic least-squares model

\[
\psi = a+bX+cY+dX^2+eXY+fY^2.
\]

The stationary point is obtained from the fitted gradient. Hyperbolic points satisfy

\[
\det H < 0,
\]

where

\[
H =
\begin{pmatrix}
2d & e\\
e & 2f
\end{pmatrix}.
\]

## 3. Saddle levels

Each saddle defines a candidate separatrix level

\[
\psi_s = \psi(x_s,y_s).
\]

Saddles whose levels agree within `level_tol` are processed together. This allows a single
heteroclinic network to contain more than one saddle.

## 4. Contour reconstruction

The level set \(\psi=\psi_s\) is extracted with marching squares. Ambiguous four-crossing cells
use a bilinear/asymptotic-style decider.

## 5. Local saddle topology

A small disk is removed around each saddle. The local quadratic form

\[
H_{11}X^2+2H_{12}XY+H_{22}Y^2=0
\]

provides the four asymptotic separatrix directions. Intersections of the contour with the disk
are matched to these directions and connected back to the saddle.

## 6. Graph assembly

Contour segments are converted to an undirected graph. Connected components containing at least
one saddle become separatrix networks. Paths between special vertices are stored as branches.

A branch can be:

- boundary-going;
- homoclinic (`saddle_start == saddle_end > 0`);
- heteroclinic (`saddle_start != saddle_end`, both nonzero).

The algorithm does not assume a fixed number of saddles, networks, or branches.


## Rectilinear nonuniform grids

The x and y coordinates may be nonuniform provided each axis is strictly increasing. Interior first derivatives
use the three-point second-order formula for unequal adjacent spacings. Marching squares uses the actual cell
coordinates, so contour interpolation remains cell-local on a nonuniform rectilinear grid.

## Graph construction complexity

Contour endpoints are quantized into tolerance-sized spatial hash cells. Each endpoint only searches its own
cell and the eight neighbours. This replaces the previous full scan over all existing vertices and makes graph
assembly approximately linear for ordinary contour data.
