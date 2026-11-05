#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from run_gmm import run_gmm

p=argparse.ArgumentParser()
p.add_argument('--output', required=True, type=Path)
p.add_argument('--prefix', default='')
p.add_argument('--locus', required=True)
a=p.parse_args()
m=re.fullmatch(r'(?:chr)?(\d+):(\d+)-(\d+)', a.locus, re.I)
if not m: p.error('Invalid locus')
chrom, start, end=m.groups()
name=f'chr{chrom}_{start}_{end}'
if a.prefix: name=f'{a.prefix}_{name}'
out=a.output.expanduser().resolve()/name
run_gmm(out/'merged_distances.tsv', out/'reference_interval.tsv', out/'gmm')
