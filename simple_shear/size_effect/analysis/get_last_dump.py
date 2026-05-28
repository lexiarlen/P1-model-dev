#!/usr/bin/env python3
# gets last dump from lammps data for making figure 3a
"""
Extract the last LAMMPS dump block from every .lammps file inside each LXXX folder,
preserve the seed in the filename, and rename 'v_r' to 'r' in the final dump header.

Expected layout:
parent/
    L50/
        *.lammps
    L60/
        *.lammps
    L100/
        *.lammps
    ...

Creates:
parent/
    last_dumps/
        L50__originalfilename.lammps
        L60__originalfilename.lammps
        ...

Only the last dump block is written.
If the header line contains 'v_r', it is replaced with 'r' for easier visualization with the OVITO GUI.
"""

from pathlib import os, Path

MARKER = b"ITEM: TIMESTEP"
ATOMS_HEADER_PREFIX = b"ITEM: ATOMS "
CHUNK_SIZE = 8 * 1024 * 1024  # 8 MB


def find_last_marker_offset(file_path: Path, marker: bytes = MARKER):
    """Find the byte offset of the last valid 'ITEM: TIMESTEP' marker."""
    file_size = file_path.stat().st_size
    if file_size == 0:
        return None

    overlap = len(marker) - 1

    with file_path.open("rb") as f:
        pos = file_size
        tail = b""

        while pos > 0:
            read_size = min(CHUNK_SIZE, pos)
            pos -= read_size
            f.seek(pos)
            chunk = f.read(read_size)

            data = chunk + tail
            idx = data.rfind(marker)

            while idx != -1:
                absolute = pos + idx

                if absolute == 0:
                    return absolute

                f.seek(absolute - 1)
                prev = f.read(1)
                if prev in (b"\n", b"\r"):
                    return absolute

                idx = data.rfind(marker, 0, idx)

            tail = data[:overlap]

    return None


def replace_vr_in_atoms_header(line: bytes) -> bytes:
    """
    Replace token 'v_r' with 'r' only in the ITEM: ATOMS header line.
    Leaves all data lines unchanged.
    """
    if line.startswith(ATOMS_HEADER_PREFIX):
        parts = line.rstrip(b"\r\n").split()
        parts = [b"radius" if p == b"v_r" else p for p in parts]

        if line.endswith(b"\r\n"):
            return b" ".join(parts) + b"\r\n"
        elif line.endswith(b"\n"):
            return b" ".join(parts) + b"\n"
        else:
            return b" ".join(parts)

    return line


def copy_from_offset_with_header_edit(src: Path, dst: Path, offset: int):
    """
    Copy src[offset:] to dst.
    If an 'ITEM: ATOMS ...' header line is encountered, replace 'v_r' -> 'r'.
    """
    with src.open("rb") as fin, dst.open("wb") as fout:
        fin.seek(offset)

        while True:
            line = fin.readline()
            if not line:
                break

            if line.startswith(ATOMS_HEADER_PREFIX):
                line = replace_vr_in_atoms_header(line)

            fout.write(line)


def process_all(parent_dir="."):
    parent = Path(parent_dir).resolve()
    output_dir = parent / "last_dumps"
    output_dir.mkdir(exist_ok=True)

    for subfolder in sorted(parent.iterdir()):
        if not subfolder.is_dir():
            continue
        if subfolder.name == "last_dumps":
            continue
        if not subfolder.name.startswith("L"):
            continue

        print(f"\nChecking folder: {subfolder.name}")

        for file_path in sorted(subfolder.glob("*.lammps")):
            try:
                print(f"  Processing: {file_path.name}")
                offset = find_last_marker_offset(file_path)

                if offset is None:
                    print("    Skipped: no 'ITEM: TIMESTEP' found")
                    continue

                out_name = f"{subfolder.name}__{file_path.name}"
                out_path = output_dir / out_name

                copy_from_offset_with_header_edit(file_path, out_path, offset)
                print(f"    Saved: {out_path.name}")

            except Exception as e:
                print(f"    Error processing {file_path.name}: {e}")

if __name__ == "__main__":
    process_all(os.environ.get("DATA_DIR", "data/simple_shear/size_effect/ensemble10_smno/"))