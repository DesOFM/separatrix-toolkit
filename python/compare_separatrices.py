#!/usr/bin/env python3
import argparse
from separatrix_viz.compare import compare_files

p = argparse.ArgumentParser()
p.add_argument('reference')
p.add_argument('candidate')
p.add_argument('--reference-schema', choices=['xy','toolkit'], default='xy')
p.add_argument('--candidate-schema', choices=['xy','toolkit'], default='toolkit')
a = p.parse_args()
r = compare_files(a.reference, a.candidate, reference_schema=a.reference_schema, candidate_schema=a.candidate_schema)
for k, v in r.items():
    print(f'{k}: {v:.8e}')
