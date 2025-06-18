#!/usr/bin/env python3
"""
Wyckoff Positions Fortran DB Generator
-------------------------------------
• Needs: numpy, spglib
• Outputs: Wyckoff positions formatted for Fortran arrays
"""
import pathlib
import numpy as np
import spglib

N_GROUPS = 530
ROT_MAX = 192

lines = []
n_sites_lines = []

for hall_number in range(1, N_GROUPS + 1):
    symmetry = spglib.get_symmetry_from_database(hall_number)

    if symmetry is None or 'rotations' not in symmetry or 'translations' not in symmetry:
        print(f"Warning: Missing symmetry data for Hall number {hall_number}")
        continue  # Skip invalid entries

    rotations = symmetry['rotations']
    translations = symmetry['translations']
    n_ops = len(rotations)

    n_sites_lines.append(f"n_sites({hall_number}) = {n_ops}")

    for op_index, (R, t) in enumerate(zip(rotations, translations), 1):
        t12 = np.rint(t * 12).astype(int)
        for row in range(3):
            lines.append(
                f"wyckoff_sites({row+1},1:4,{hall_number},{op_index}) = [{R[row,0]}, {R[row,1]}, {R[row,2]}, {t12[row]}]"
            )

# Output the Fortran-compatible file
out = pathlib.Path("wyckoff_positions_fortran.f90")
out.write_text("\n".join(lines + ["", "! Number of sites per Hall symbol"] + n_sites_lines))
print(f"✓ Wrote {out} ({len(lines)} lines + {len(n_sites_lines)} n_sites entries)")
