# Hydrodynamic validation model

The separatrix-reconstruction algorithm implemented in this repository was developed independently.

The motivating hydrodynamic validation case uses a three-layer vortex-flow model from the
contour-dynamics framework described by M. A. Sokolovskiy. The hydrodynamic model is used only
to generate the streamfunction field supplied to the separatrix algorithm; it is not part of the
reusable separatrix-reconstruction method.

The two model references used for this validation context are:

1. Соколовский М. А. Моделирование трехслойных вихревых движений в океане методом контурной
   динамики // Известия АН СССР. Физика атмосферы и океана. 1991. Т. 27, № 5. С. 550–562.

2. Sokolovskiy M. A. Stability Analysis of the Axisymmetric Three-Layered Vortex Using Contour
   Dynamics Method // Computational Fluid Dynamics Journal. 1997. Vol. 6, No. 2. P. 133–156.

Machine-readable BibTeX records are provided in [`REFERENCES.bib`](../REFERENCES.bib).

## Separation of responsibilities

Conceptually, the workflow is

```text
Sokolovskiy three-layer model
            |
            v
      streamfunction field
          psi(x,y)
            |
            v
  Separatrix Toolkit algorithm
            |
            v
 saddle points + separatrix networks
```

The citation of the hydrodynamic model therefore documents the provenance of one validation field.
It does not imply that the separatrix-reconstruction algorithm or its implementation was taken from
those publications.
