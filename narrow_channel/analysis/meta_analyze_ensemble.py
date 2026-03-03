import numpy as np
from pathlib import Path
import os
import argparse
import xarray as xr

# --------------------
# global vars
# --------------------
#vari = "gamma"

# parameter indices in THETA
KB_IDX  = 0
EPS_IDX = 1

parser = argparse.ArgumentParser()
parser.add_argument("--J", type=int, required=True)
parser.add_argument("--parent_dir_name", type=str, required=False)
parser.add_argument("--power", type=float, default=1)
parser.add_argument("--variable", type=str, required=True)
parser.add_argument("--natoms", type=int, required=False)
args = parser.parse_args()
J = args.J
vari = args.variable
POWER = 1/float(args.power)
parent_dir_name = args.parent_dir_name if args.parent_dir_name is not None else "ensemble"
natoms = args.natoms if args.natoms is not None else 3000

output_dir = Path(f"/scratch/groups/earlew/arlenlex/P1-EKI-nares/EKI/{parent_dir_name}/")

if vari == "u":
    varis = np.linspace(5, 30, 11)
if vari == "rho":
    varis = np.linspace(600, 1600, 11)
if vari == "rho_a":
    varis = np.linspace(0.7, 1.7, 11)
if vari == "rho_o":
    varis = np.linspace(800, 1300, 11)
if vari == "C_a":
    varis = np.linspace(0.0005, 0.003, 11)
if vari == "T_max":
    varis = np.linspace(20, 40, 11)
if vari == "gamma":
    varis = np.linspace(1e8, 1e10, 11)
if vari == "damping_ratio":
    varis = np.linspace(0.5, 4.5, 11)
if vari == "C":
    varis = np.linspace(30e3,75e3,10)
if vari == "domain_scaling":
    varis = np.linspace(0.6, 2.0, 11)
# --------------------
# io helpers
# --------------------
def load_outputs(state_dir: Path):
    """Load THETA (2 x J), Y (2 x J), SIG (2 x J) from one directory."""
    ftheta = os.path.join(state_dir, "theta.tsv")
    THETA = np.array(np.loadtxt(ftheta), dtype=float).T           # shape (2, J)
    print(f'theta shape = {THETA.shape}')

    fy = os.path.join(state_dir, "Y.npy")
    Y = np.load(fy)                                               # shape (2, J)
    print(f'y shape = {Y.shape}')

    fsig = os.path.join(state_dir, "SIG.npy")
    SIG = np.load(fsig)                                           # shape (2, J)
    print(f'sigma shape = {SIG.shape}')

    return THETA, Y, SIG

def mask_arrays(THETA, Y, SIG):
    """Mask arrays by NaN where |Y|==999 and drop columns with any NaN in Y."""
    Y_masked = Y.astype(float)
    Y_masked[np.abs(Y_masked) == 999] = np.nan

    THETA_masked = THETA.astype(float).copy()
    SIG_masked   = SIG.astype(float).copy()

    colmask = np.isnan(Y_masked).any(axis=0)  # columns to drop

    Y_masked[:, colmask]     = np.nan
    SIG_masked[:, colmask]   = np.nan
    THETA_masked[:, colmask] = np.nan
    return THETA_masked, Y_masked, SIG_masked

# --------------------
# algebraic fit (no intercept)
# --------------------
def algebraic_fit_both_peaks(Z_masked, THETA_masked, power=0.5):
    """
    For each peak (0,1), fit z = A * (k_b*eps)^power (no intercept).
    Solve A by linear least squares on x = (k_b*eps)^power.
    Returns arrays (2,) for A, s (=power), and r2 (linear space).
    """
    kb  = THETA_masked[KB_IDX, :].ravel()
    eps = THETA_masked[EPS_IDX, :].ravel()
    p   = kb * eps  # predictor

    As  = np.full(2, np.nan)
    ss  = np.full(2, np.nan)   # store input power
    r2s = np.full(2, np.nan)

    for peak_idx in (0, 1):
        z = Z_masked[peak_idx, :].ravel()

        # finite values; require p >= 0 for fractional powers
        valid = np.isfinite(p) & np.isfinite(z) & (p >= 0)
        if not np.any(valid):
            continue

        pv, zv = p[valid], z[valid]
        x = pv ** power

        denom = np.sum(x * x)
        if denom <= 0:
            continue

        A = np.sum(x * zv) / denom
        zhat = A * x
        ss_res = np.sum((zv - zhat) ** 2)
        ss_tot = np.sum((zv - np.mean(zv)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

        As[peak_idx]  = A
        ss[peak_idx]  = power
        r2s[peak_idx] = r2

    return As, ss, r2s

# --------------------
# regime delineation (DP segmentation on s = log10(kb*eps))
# --------------------
def delineate_regime_boundaries(THETA, Y, K=4):
    """
    Build regime labels from Y (==999 encodes missing peak),
    sort by s = log10(kb*eps), and find K segments that minimize
    misclassification errors. Return:
      cuts_ke: (K-1,) thresholds in kb*eps
      acc:     float accuracy of the K-seg labeling
      counts:  (4,) counts per regime [0:no frac, 1:stable, 2:no arch, 3:unstable]
    """
    y1 = Y[0, :].astype(float)
    y2 = Y[1, :].astype(float)
    m1 = (y1 == 999)
    m2 = (y2 == 999)

    kb  = THETA[KB_IDX, :].astype(float).ravel()
    eps = THETA[EPS_IDX, :].astype(float).ravel()

    valid = np.isfinite(kb) & np.isfinite(eps) & (kb > 0) & (eps > 0)
    if not np.any(valid) or valid.sum() < K:
        return np.full(K-1, np.nan), np.nan, np.zeros(4, dtype=int)

    kb, eps = kb[valid], eps[valid]
    m1, m2  = m1[valid], m2[valid]

    # regimes
    regime = np.full(kb.shape, 3, dtype=int)
    regime[m1 &  m2] = 0
    regime[~m1 & m2] = 1
    regime[m1 & ~m2] = 2

    # sort by s = log10(kb*eps)
    s = np.log10(kb * eps).astype(float)
    idx = np.argsort(s)
    s, y = s[idx], regime[idx]
    n = len(s)

    # encode labels to 0..C-1
    labels, y_enc = np.unique(y, return_inverse=True)
    C = len(labels)

    # prefix counts
    P = np.zeros((C, n+1), dtype=int)
    for c in range(C):
        P[c, 1:] = np.cumsum((y_enc == c).astype(int))

    def seg_cost(i, j):
        counts = P[:, j] - P[:, i]
        seglen = (j - i)
        return seglen - counts.max()

    # DP
    dp = np.full((K+1, n+1), np.inf)
    bp = np.full((K+1, n+1), -1, dtype=int)
    dp[0, 0] = 0.0

    # precompute interval costs
    cost = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n+1):
            cost[i, j-1] = seg_cost(i, j)

    for k in range(1, K+1):
        for j in range(1, n+1):
            best_i = -1
            best_v = np.inf
            for i in range(k-1, j):  # at least 1 point per segment
                v = dp[k-1, i] + cost[i, j-1]
                if v < best_v:
                    best_v, best_i = v, i
            dp[k, j] = best_v
            bp[k, j] = best_i

    # recover cut indices
    cuts_idx = []
    k, j = K, n
    while k > 0:
        i = bp[k, j]
        cuts_idx.append(i)
        j, k = i, k-1
    cuts_idx = sorted(cuts_idx)[1:]  # internal boundaries (K-1)

    # midpoints in s, then back to kb*eps
    cuts_s = []
    for i in cuts_idx:
        if i == 0 or i >= n:
            cuts_s.append(np.nan)
        else:
            cuts_s.append(0.5 * (s[i-1] + s[i]))
    cuts_ke = 10.0 ** np.array(cuts_s)

    # compute accuracy
    pred = np.empty(n, dtype=int)
    edges = [0] + cuts_idx + [n]
    for a, b in zip(edges[:-1], edges[1:]):
        cnts = P[:, b] - P[:, a]
        maj = labels[np.argmax(cnts)]
        pred[a:b] = maj
    acc = (pred == y).mean()

    counts = np.array([(regime == r).sum() for r in (0, 1, 2, 3)], dtype=int)
    return cuts_ke, acc, counts

# --------------------
# collect across vari
# --------------------
all_A, all_s, all_r2s = [], [], []

all_sigmin, all_sigmax, all_N = [], [], []
vals_used = []
all_ke_bounds = []
all_acc = []
all_counts = []

unit = "kg/m^3" if vari == "rho" else "m/s"

if vari == "gamma":
    string_varis = ['100e6', '109e7', '208e7', '307e7', '406e7', '505e7', '604e7',
                 '703e7', '802e7', '901e7', '100e8']
if vari == "C":
    string_varis = ['30e3', '35e3', '40e3','45e3', '50e3', '55e3', '60e3', '65e3', '70e3', '75e3']


for i, v in enumerate(varis):
    state_dir = output_dir / f"N{natoms}_{vari}{int(v)}" / "state"
    if vari == 'gamma':
        state_dir = output_dir / f"N{natoms}_{vari}{string_varis[i]}" / "state"
    if vari == "rho_o":
        state_dir = output_dir / f"N{natoms}_rho{int(v)}" / "state"
    if vari == "rho_a":
        state_dir = output_dir / f"N{natoms}_rho{int(v*10)}" / "state"
    if vari == "damping_ratio":
        state_dir = output_dir / f"N{natoms}_damping_ratio{int(v*10)}" / "state"
    if vari == "C_a":
        state_dir = output_dir / f"N{natoms}_C_a{int(v*100000)}" / "state"
    if vari == "domain_scaling":
        state_dir = output_dir / f"N{natoms}_domain_scaling{np.round(v*100):.0f}" / "state"
    if vari == "C":
        state_dir = output_dir / f"N{natoms}_{vari}{string_varis[i]}" / "state"
    
    if not state_dir.is_dir():
        continue

    print(f'working on {vari} = {v} {unit}')
    THETA, Y, SIG = load_outputs(state_dir)
    THETA_masked, Y_masked, SIG_masked = mask_arrays(THETA, Y, SIG)

    # N = number of columns where both kb and eps are finite (after masking)
    N_theta = int(np.isfinite(THETA_masked).all(axis=0).sum())
    all_N.append(N_theta)

    # σ = kb*eps spread (pre-mask but finite)
    sig = THETA_masked[KB_IDX, :] * THETA_masked[EPS_IDX, :]
    valid_sig = np.isfinite(sig)
    all_sigmin.append(np.nanmin(sig[valid_sig]) if np.any(valid_sig) else np.nan)
    all_sigmax.append(np.nanmax(sig[valid_sig]) if np.any(valid_sig) else np.nan)

    # Regime boundaries from raw Y/THETA
    ke_bounds, acc, counts = delineate_regime_boundaries(THETA, Y, K=4)
    all_ke_bounds.append(ke_bounds)
    all_acc.append(acc)
    all_counts.append(counts)

    # Algebraic fits (only)
    A, s, r2s = algebraic_fit_both_peaks(SIG_masked, THETA_masked, power=POWER)
    all_A.append(A)
    all_s.append(s)          # store the chosen power per peak
    all_r2s.append(r2s)

    vals_used.append(float(v))

# --------------------
# pack to xarray
# --------------------
if len(vals_used) > 0:
    peak = np.array([1, 2])
    boundary = np.arange(1, 4)
    regime = np.array([0, 1, 2, 3])

    all_A       = np.vstack(all_A)
    all_s       = np.vstack(all_s)
    all_r2s     = np.vstack(all_r2s)

    all_sigmin    = np.asarray(all_sigmin)
    all_sigmax    = np.asarray(all_sigmax)
    all_N         = np.asarray(all_N, dtype=np.int32)
    all_ke_bounds = np.asarray(all_ke_bounds, dtype=float)
    all_acc       = np.asarray(all_acc, dtype=float)
    all_counts    = np.asarray(all_counts, dtype=np.int32)

    VAR = np.asarray(vals_used)

    ds = xr.Dataset(
        data_vars=dict(
            # Algebraic fit outputs
            A_algebraic=([vari, "peak"], all_A),          # fitted A (per peak)
            power_used=([vari, "peak"], all_s),           # the input 'power' duplicated per peak
            r2_algebraic=([vari, "peak"], all_r2s),

            # Spread and counts
            sig_min=([vari], all_sigmin),
            sig_max=([vari], all_sigmax),
            N=([vari], all_N),

            # Regime delineation on kb*eps
            ke_bounds=([vari, "boundary"], all_ke_bounds),
            regime_acc=([vari], all_acc),
            regime_counts=([vari, "regime"], all_counts),
        ),
        coords={
            vari: (vari, VAR),
            "peak": ("peak", peak),
            "boundary": ("boundary", boundary),
            "regime": ("regime", regime),
        },
        attrs=dict(
            description="Algebraic fits z = A*(k_b*eps)^power for both peaks, plus regime thresholds in k_b*eps.",
            algebraic_model=f"z = A_algebraic * (k_b*eps)^(power_used); power_used == {POWER} (user input). Intercept fixed to 0.",
            regimes="0=no fracture (999,999); 1=stable arch (x,999); 2=no arch (999,x); 3=unstable arch (x,x)",
            segmentation="Dynamic programming on sorted s=log10(k_b*eps) to find 4 segments (3 boundaries) minimizing misclassification.",
            note="Masked |Y|==999 for fits only. Regime delineation uses raw Y==999 flags.",
        ),
    )

    out_nc = output_dir / "algebraic_fits.nc"
    ds.to_netcdf(out_nc)
    print(f"Saved: {out_nc}")
else:
    print(f"No matching state directories found in {output_dir} Nothing saved.")
