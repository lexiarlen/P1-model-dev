#!/usr/bin/env python3

import os
import re
import json
import argparse


def hms_to_seconds(hms):
    """
    Convert HH:MM:SS string to seconds.
    """
    h, m, s = map(int, hms.split(":"))
    return h * 3600 + m * 60 + s


def parse_lammps_runs_to_dict(text):
    """
    Parse LAMMPS output text and return a dictionary grouped by atom count.

    Output format:
    {
        n_atoms: {
            "total_wall_time_s": [...],
            "n_cpus": [...]
        }
    }
    """

    created_atoms_re = re.compile(r"Created\s+(\d+)\s+atoms")

    cpu_re = re.compile(
        r"CPU use with\s+(\d+)\s+MPI tasks\s+x\s+(\d+)\s+OpenMP threads"
    )

    loop_re = re.compile(
        r"Loop time of\s+[0-9.eE+-]+\s+on\s+(\d+)\s+procs\s+for\s+\d+\s+steps\s+with\s+(\d+)\s+atoms"
    )

    wall_re = re.compile(
        r"Total wall time:\s*([0-9]+:[0-9]{2}:[0-9]{2})"
    )

    runs_by_atoms = {}

    current_atoms = None
    current_cpus = None

    for line in text.splitlines():
        line = line.strip()

        m_atoms = created_atoms_re.search(line)
        if m_atoms:
            current_atoms = int(m_atoms.group(1))
            current_cpus = None
            continue

        m_cpu = cpu_re.search(line)
        if m_cpu:
            n_mpi = int(m_cpu.group(1))
            n_omp = int(m_cpu.group(2))
            current_cpus = n_mpi * n_omp
            continue

        m_loop = loop_re.search(line)
        if m_loop and current_cpus is None:
            current_cpus = int(m_loop.group(1))

            if current_atoms is None:
                current_atoms = int(m_loop.group(2))

            continue

        m_wall = wall_re.search(line)
        if m_wall and current_atoms is not None:
            if current_atoms not in runs_by_atoms:
                runs_by_atoms[current_atoms] = {
                    "total_wall_time_s": [],
                    "n_cpus": [],
                }

            runs_by_atoms[current_atoms]["total_wall_time_s"].append(
                hms_to_seconds(m_wall.group(1))
            )

            runs_by_atoms[current_atoms]["n_cpus"].append(current_cpus)

    return runs_by_atoms


def merge_runs_dicts(dict_list):
    """
    Merge several runs dictionaries into one.
    """

    merged = {}

    for d in dict_list:
        for n_atoms, vals in d.items():
            if n_atoms not in merged:
                merged[n_atoms] = {
                    "total_wall_time_s": [],
                    "n_cpus": [],
                }

            merged[n_atoms]["total_wall_time_s"].extend(vals["total_wall_time_s"])
            merged[n_atoms]["n_cpus"].extend(vals["n_cpus"])

    return merged


def save_runs_dict(runs_dict, out_path):
    """
    Save runs_dict as JSON with plain Python int keys/values.
    """

    runs_dict_clean = {
        int(n_atoms): {
            "total_wall_time_s": [int(x) for x in vals["total_wall_time_s"]],
            "n_cpus": [
                int(x) if x is not None else None
                for x in vals["n_cpus"]
            ],
        }
        for n_atoms, vals in runs_dict.items()
    }

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w") as f:
        json.dump(runs_dict_clean, f, indent=2)

    print(f"Saved runs_dict to: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Parse LAMMPS output files and save runs_dict.json."
    )

    parser.add_argument(
        "--base-dir",
        required=True,
        help="Directory containing LAMMPS .out files.",
    )

    parser.add_argument(
        "--files",
        nargs="+",
        required=True,
        help="LAMMPS output file names relative to base-dir.",
    )

    parser.add_argument(
        "--out-path",
        required=True,
        help="Path where runs_dict.json should be saved.",
    )

    args = parser.parse_args()

    all_runs_dicts = []

    for fname in args.files:
        log_path = os.path.join(args.base_dir, fname)

        print(f"Parsing: {log_path}")

        if not os.path.exists(log_path):
            raise FileNotFoundError(f"Could not find file: {log_path}")

        with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()

        parsed = parse_lammps_runs_to_dict(text)
        all_runs_dicts.append(parsed)

    runs_dict = merge_runs_dicts(all_runs_dicts)

    print("\nParsed runs:")
    for n_atoms in sorted(runs_dict.keys()):
        n_runs = len(runs_dict[n_atoms]["total_wall_time_s"])
        print(f"{n_atoms:>10} atoms | {n_runs:>3} runs")

    save_runs_dict(runs_dict, args.out_path)


if __name__ == "__main__":
    main()