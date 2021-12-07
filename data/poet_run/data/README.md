As an example, here is a description of the files in DFT_Cu/NVT_FCC_300K

# NVT_FCC_300K
directory with snapshots of density functional theory molecular dynamics at 300K on the NVT ensemble

# NVT_FCC_1400K
directory with snapshots of density functional theory molecular dynamics at 1400K on the NVT ensemble

# NPT_FCC_1400K
directory with snapshots of density functional theory molecular dynamics at 1400K and 100kPa on the NPT ensemble

# lowIndexSurfaces
directory with thirteen low-index surfaces of fcc copper

# 0K
directory with data for copper at zero Kelvin

# F_coord.data
cartesian coordinates (Angstroms) (first three columns, x, y, and z, respectively) and forces (eV/Angstrom) (next three columns, x, y, and z, respectively) of each site in the structure. For example, the first thirty-two rows correspond to the first structure (i.e., 1.vasp), the second thirty-two correspond to the second structure (i.e., 2.vasp)

# F_Cell.data
negative components of the virial stress tensor in units of kB. There is one row for each structure in the directory. Each column corresponds to the components XX, YY, ZZ, XY, YZ, and ZX, respectively

# E0.data
energies corresponding the each structure in the directory. For example, the first row corresponds to the first structure (i.e., 1.vasp), the second row corresponds to the second structure (i.e., 2.vasp)

# 1.vasp
structure in POSCAR format (for example, 1.vasp corresponds to the first snapshot in the case of molecular dynamics data)

# TargetProperties_0K
file with data for copper at zero Kelvin