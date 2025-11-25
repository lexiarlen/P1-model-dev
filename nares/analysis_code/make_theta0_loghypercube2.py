#!/usr/bin/env python3
"""

Writes a TSV with N rows and columns: kb kc eps
Each parameter is sampled from a log hypercube with user-specified min and max.

Example:
  python make_theta0_loghypercube.py \
    --kb-min 1e2  --kb-max 1e4 \
    --eps-min 1e-3 --eps-max 1e-1 \
    --J 10 --seed 42 --out state/theta_iter0.tsv
"""

import argparse
from pathlib import Path
import numpy as np
from scipy.stats import qmc

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--kb-min", type=float, required=True)
    p.add_argument("--kb-max", type=float, required=True)
    p.add_argument("--eps-min", type=float, required=True)
    p.add_argument("--eps-max", type=float, required=True)
    p.add_argument("--J", type=int, default=10, help="ensemble size")
    p.add_argument("--seed", type=int, default=428)
    p.add_argument("--out", type=str, default="state/theta_iter0.tsv")
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    Path(Path(args.out).parent).mkdir(parents=True, exist_ok=True)

    sampler = qmc.LatinHypercube(d=2, rng=rng)
    sample = sampler.random(n=args.J)
    upper_bounds = [np.log10(args.kb_max), np.log10(args.eps_max)]
    lower_bounds = [np.log10(args.kb_min), np.log10(args.eps_min)]
    theta = 10**(qmc.scale(sample, lower_bounds, upper_bounds))

    np.savetxt(args.out, theta, fmt="%.8g", header="kb\teps")
    print(f"Created initial theta and saved it to {args.out}")

if __name__ == "__main__":
    main()
