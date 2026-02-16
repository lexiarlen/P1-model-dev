# convert numpy array to netcdf for data intepretability
#!/usr/bin/env python3
"""
Build GE_results.nc from saved .npy arrays.

Modes (auto-detected):
1) RATES mode: gamma_dot.npy present (kb, eps constants; arrays per shear-rate)
2) CASES mode: one result per kb/eps (or damping_ratio) case

Now also ingests damping_ratio.npy:
  - rates mode: length-nrate -> coord "damping_ratio" on dim "rate";
                0-D/scalar   -> 0-D DataArray "damping_ratio"
  - cases mode: length-ncase -> data_var "damping_ratio" on dim "case" (can be NaN)
"""

import os
import argparse
import numpy as np
import xarray as xr
from datetime import datetime

def _maybe(base, name):
    p = os.path.join(base, name)
    return np.load(p) if os.path.exists(p) else None

def _as_scalar_da(name, value, attrs=None):
    da = xr.DataArray(value)
    da.name = name
    if attrs: da.attrs.update(attrs)
    return da

def _is_scalar(arr):
    a = np.asarray(arr)
    return a.shape == () or a.ndim == 0 or (a.ndim == 1 and a.size == 1)

def build_rates_mode(base, args):
    # required
    kb   = np.load(os.path.join(base, "kb.npy"))
    eps  = np.load(os.path.join(base, "eps.npy"))
    rates = np.load(os.path.join(base, "gamma_dot.npy"))
    Gm   = np.load(os.path.join(base, "G_mean.npy"))
    Em   = np.load(os.path.join(base, "E_mean.npy"))

    # optional
    Gs = _maybe(base, "G_std.npy")
    Es = _maybe(base, "E_std.npy")
    pg = _maybe(base, "peak_gamma_mean.npy")
    pt = _maybe(base, "peak_tau_mean.npy")
    nf = _maybe(base, "n_files.npy")
    dr = _maybe(base, "damping_ratio.npy")

    nrate = rates.shape[0]
    for arr, nm in [(Gm,"G_mean"), (Em,"E_mean"), (Gs,"G_std"), (Es,"E_std"),
                    (pg,"peak_gamma_mean"), (pt,"peak_tau_mean"), (nf,"n_files")]:
        if arr is None: continue
        if arr.shape[0] != nrate:
            raise ValueError(f"{nm} length {arr.shape[0]} != gamma_dot length {nrate}")

    ds = xr.Dataset(
        data_vars={
            "G_mean": ("rate", Gm),
            "E_mean": ("rate", Em),
        },
        coords={
            "rate": ("rate", rates, {"units": "1/s", "long_name": "shear rate (gamma_dot)"}),
        },
        attrs=dict(
            created=datetime.utcnow().isoformat() + "Z",
            source="make_nc_from_numpy (rates mode)",
            nu=args.nu, phi=args.phi,
            L0=args.L0, run_time=args.run_time, num_dumps=args.num_dumps,
            notes="kb, eps constants across rates; per-rate moduli.",
        ),
    )

    if Gs is not None: ds["G_std"] = ("rate", Gs)
    if Es is not None: ds["E_std"] = ("rate", Es)
    if pg is not None: ds["peak_gamma_mean"] = ("rate", pg)
    if pt is not None: ds["peak_tau_mean"] = ("rate", pt)
    if nf is not None: ds["n_files"] = ("rate", nf.astype(int))

    # kb, eps as true scalars (0-D)
    ds["kb"]  = _as_scalar_da("kb", float(np.asarray(kb)),  {"units": "N/m"})
    ds["eps"] = _as_scalar_da("eps", float(np.asarray(eps)), {"units": "-"})

    # damping ratio handling
    if dr is not None:
        if _is_scalar(dr):
            ds["damping_ratio"] = _as_scalar_da("damping_ratio", float(np.asarray(dr)))
        else:
            if dr.shape[0] != nrate:
                raise ValueError(f"damping_ratio length {dr.shape[0]} != gamma_dot length {nrate}")
            ds = ds.assign_coords(damping_ratio=("rate", dr))

    # convenience: sorted by rate
    ds = ds.sortby("rate")
    return ds

def build_cases_mode(base, args):
    # required
    kb   = np.load(os.path.join(base, "kb.npy"))
    eps  = np.load(os.path.join(base, "eps.npy"))
    prod = np.load(os.path.join(base, "product_kb_eps.npy"))
    Gm   = np.load(os.path.join(base, "G_mean.npy"))
    Em   = np.load(os.path.join(base, "E_mean.npy"))

    # optional
    Gs = _maybe(base, "G_std.npy")
    Es = _maybe(base, "E_std.npy")
    pg = _maybe(base, "peak_gamma_mean.npy")
    pt = _maybe(base, "peak_tau_mean.npy")
    nf = _maybe(base, "n_files.npy")
    dr = _maybe(base, "damping_ratio.npy")   # NEW

    ncase = kb.shape[0]
    for arr, nm in [(eps,"eps"), (prod,"product_kb_eps"), (Gm,"G_mean"), (Em,"E_mean"),
                    (Gs,"G_std"), (Es,"E_std"), (pg,"peak_gamma_mean"),
                    (pt,"peak_tau_mean"), (nf,"n_files")]:
        if arr is None: continue
        if arr.shape[0] != ncase:
            raise ValueError(f"{nm} length {arr.shape[0]} != kb length {ncase}")
    if dr is not None and not _is_scalar(dr) and dr.shape[0] != ncase:
        raise ValueError(f"damping_ratio length {dr.shape[0]} != kb length {ncase}")

    ds = xr.Dataset(
        data_vars={
            "kb": ("case", kb, {"units": "N/m"}),
            "eps": ("case", eps, {"units": "-"}),
            "product_kb_eps": ("case", prod, {"units": "N/m"}),
            "G_mean": ("case", Gm),
            "E_mean": ("case", Em),
        },
        coords={"case": np.arange(ncase, dtype=int)},
        attrs=dict(
            created=datetime.utcnow().isoformat() + "Z",
            source="make_nc_from_numpy (cases mode)",
            nu=args.nu, phi=args.phi, shear_rate=args.shear_rate,
            L0=args.L0, run_time=args.run_time, num_dumps=args.num_dumps,
            notes="One result per case. damping_ratio present if provided.",
        ),
    )

    if Gs is not None: ds["G_std"] = ("case", Gs)
    if Es is not None: ds["E_std"] = ("case", Es)
    if pg is not None: ds["peak_gamma_mean"] = ("case", pg)
    if pt is not None: ds["peak_tau_mean"] = ("case", pt)
    if nf is not None: ds["n_files"] = ("case", nf.astype(int))
    if dr is not None:
        if _is_scalar(dr):
            # Treat scalar as 0-D: applies to all cases (rare; mostly for completeness)
            ds["damping_ratio_scalar"] = _as_scalar_da("damping_ratio_scalar", float(np.asarray(dr)))
        else:
            ds["damping_ratio"] = ("case", dr)

    # legacy convenience coords (time, gamma) for a single shear-rate path
    t = np.linspace(0.0, args.run_time, args.num_dumps + 1)
    gamma = args.shear_rate * t
    ds["time"] = ("time", t); ds["time"].attrs["units"] = "s"
    ds["gamma"] = ("time", gamma); ds["gamma"].attrs["description"] = "shear strain = shear_rate * time"

    # Sorting: keep legacy by kb; if kb identical and damping_ratio exists, also sort by it.
    ds = ds.sortby("kb")
    if "damping_ratio" in ds and float(np.nanstd(ds["kb"].values)) == 0.0:
        ds = ds.sortby("damping_ratio")

    return ds

def main():
    ap = argparse.ArgumentParser(description="Build NetCDF from saved .npy arrays (includes damping_ratio if present)")
    ap.add_argument("--base-dir", required=True, help="Folder containing the .npy files")
    # optional reference values to store in attrs
    ap.add_argument("--nu", type=float, default=0.3)
    ap.add_argument("--phi", type=float, default=0.71)
    ap.add_argument("--shear-rate", type=float, default=1e-7)  # used in cases mode for gamma(t)
    ap.add_argument("--L0", type=float, default=100e3)
    ap.add_argument("--run-time", type=float, default=3600.0)
    ap.add_argument("--num-dumps", type=int, default=500)
    args = ap.parse_args()

    base = os.path.abspath(args.base_dir)

    # Mode select
    if os.path.exists(os.path.join(base, "gamma_dot.npy")):
        ds = build_rates_mode(base, args)
    else:
        ds = build_cases_mode(base, args)

    out_nc = os.path.join(base, "GE_results.nc")
    ds.to_netcdf(out_nc)
    print(f"[OK] wrote {out_nc}")

if __name__ == "__main__":
    main()
