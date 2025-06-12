#!/usr/bin/env python3
"""
makeSpaceGroup.py – create spacegroup_db.f90
--------------------------------------------
• Needs:  numpy, spglib               (pip wheels – no compiler)
• Builds: spacegroup_db.f90  (~950 kB)  with
      rotations(3,3,192,530)
      translations(3,192,530)
      n_ops(530)
      hall(530)
Run once, then   USE spacegroup_db   from Fortran.
"""
import pathlib, textwrap, numpy as np, spglib

N_GROUPS  = 530          # Hall numbers 1…530
ROT_MAX   = 192          # max ops in any space group (F-centred cubic)

HALL_FMT = '        "{:<17s}"'       # 8-space indent, fixed width 17

hall_syms, rot_lines, tran_lines, nops_lines = [], [], [], []

# ---------------------------------------------------------------
# Harvest data from spglib
# ---------------------------------------------------------------
for g in range(1, N_GROUPS + 1):
    sg   = spglib.get_spacegroup_type(g)
    ops  = spglib.get_symmetry_from_database(g)

    hall_raw  = sg.get("hall", sg.get("hall_symbol"))   # old/new key
    hall_safe = hall_raw.replace('"', '""')             # escape internal "
    hall_syms.append(HALL_FMT.format(hall_safe))

    n_ops = len(ops["rotations"])
    nops_lines.append(f"  n_ops({g}) = {n_ops}")

    for k, (R, t) in enumerate(zip(ops["rotations"], ops["translations"]), 1):
        t12 = (np.rint(t * 12)).astype(int)             # exact ints (twelfths)
        for i in range(3):
            rot_lines.append(
                f"  rotations({i+1},1:3,{k},{g}) = (/ {R[i,0]:2d}, {R[i,1]:2d}, {R[i,2]:2d} /)")
        tran_lines.append(
            f"  translations(1:3,{k},{g}) = (/ {t12[0]:2d}, {t12[1]:2d}, {t12[2]:2d} /)")

# ---------------------------------------------------------------
# Assemble the Fortran module
# ---------------------------------------------------------------
TEMPLATE = textwrap.dedent("""
    module spacegroup_db
      implicit none
      integer, parameter :: n_groups  = {ng}
      integer, parameter :: n_ops_max = {rm}

      integer, parameter :: rotations(3,3,n_ops_max,n_groups) = 0
      integer, parameter :: translations(3,n_ops_max,n_groups) = 0
      integer, parameter :: n_ops(n_groups) = 0

      character(len=17), parameter :: hall(n_groups) = (/ &
{hall_block} &
      /)

{nops_block}
{rot_block}
{tran_block}

    end module spacegroup_db
""").format(
        ng=N_GROUPS,
        rm=ROT_MAX,
        hall_block=", &\n".join(hall_syms),
        nops_block="\n".join(nops_lines),
        rot_block="\n".join(rot_lines),
        tran_block="\n".join(tran_lines)
    )

out = pathlib.Path("spacegroup_db.f90")
out.write_text(TEMPLATE)
print(f"✓ Wrote {out} ({out.stat().st_size/1024:.1f} kB)")

