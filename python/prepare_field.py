#!/usr/bin/env python3
import argparse
from separatrix_viz.io import prepare_field

p = argparse.ArgumentParser()
p.add_argument('source')
p.add_argument('output')
p.add_argument('--x-col', type=int, default=1)
p.add_argument('--y-col', type=int, default=2)
p.add_argument('--psi-col', type=int, default=3)
a = p.parse_args()
prepare_field(a.source, a.output, x_col=a.x_col, y_col=a.y_col, psi_col=a.psi_col)
