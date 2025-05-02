import os
import random
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from pymatgen.core import Lattice, Structure
from pymatgen.io.lammps.data import LammpsData
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from scipy.interpolate import interp1d

n = 4  # Supercell size

# --------------------------------------------------
# 1. Generate LSMO structure with Sr doping
# --------------------------------------------------
def generate_lsmo_structure(x=0.3, supercell=(n, n, n)):
    a = 3.9
    lattice = Lattice.cubic(a)
    species = ["La", "Mn", "O", "O", "O"]
    coords = [
        [0, 0, 0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5]
    ]
    struct = Structure(lattice, species, coords) * supercell
    la_indices = [i for i, site in enumerate(struct) if site.specie.symbol == "La"]
    num_replace = int(x * len(la_indices))
    sr_indices = random.sample(la_indices, num_replace)
    for i in sr_indices:
        struct[i] = "Sr"
    return struct

# --------------------------------------------------
# 2. Export structure to LAMMPS format
# --------------------------------------------------
def export_to_lammps(structure, filename="lsmo.lmpdat"):
    lmp_data = LammpsData.from_structure(structure, atom_style="atomic")
    lmp_data.write_file(filename)
    print(f"[INFO] LAMMPS data written to {filename}")

# --------------------------------------------------
# 3. LAMMPS input file generation
# --------------------------------------------------
def generate_lammps_input(filename="in.lsmo"):
    with open(filename, 'w') as f:
        f.write("""
units           metal
atom_style      atomic
read_data       lsmo.lmpdat

pair_style      lj/cut 2.5
pair_coeff      * * 1.0 1.0 2.5

neighbor        0.3 bin
neigh_modify    delay 5

thermo          100
min_style       cg
minimize        1e-4 1e-6 100 1000

write_data      lsmo_minimized.lmpdat
        """)
    print(f"[INFO] LAMMPS input script written to {filename}")

# --------------------------------------------------
# 4. Run LAMMPS
# --------------------------------------------------
def run_lammps(input_script="in.lsmo"):
    from lammps import lammps
    print("[INFO] Running LAMMPS...")
    lmp = lammps()
    lmp.file(input_script)
    print("[INFO] LAMMPS run complete.")

# --------------------------------------------------
# 5. Read minimized LAMMPS output
# --------------------------------------------------
def read_lammps_minimized(filename="lsmo_minimized.lmpdat"):
    if not os.path.exists(filename):
        print(f"[WARNING] {filename} not found.")
        return None
    return LammpsData.from_file(filename, atom_style="atomic").structure

# --------------------------------------------------
# 6. XRD Plot
# --------------------------------------------------
def plot_smoothed_xrd(structure, save_path="xrd_smoothed.png"):
    xrd = XRDCalculator(wavelength="CuKa")
    pattern = xrd.get_pattern(structure)
    x, y = pattern.x, pattern.y
    x_new = np.linspace(min(x), max(x), 2000)
    f = interp1d(x, y, kind='cubic', bounds_error=False, fill_value=0)
    y_smooth = f(x_new)
    y_smooth[y_smooth < 0] = 0

    plt.figure()
    plt.plot(x_new, y_smooth, color="blue", label="XRD (Smoothed)")
    plt.xlabel("2θ (degrees)")
    plt.ylabel("Intensity")
    plt.title("XRD Pattern (LAMMPS-Minimized LSMO)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"[INFO] XRD plot saved as '{save_path}'")
    plt.show()

# --------------------------------------------------
# 7. Composition Pie Chart
# --------------------------------------------------
def plot_composition(structure, save_path="composition_pie.png"):
    counts = Counter([site.specie.symbol for site in structure])
    labels = counts.keys()
    sizes = counts.values()

    plt.figure()
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140)
    plt.title("Elemental Composition (Post-Doping)")
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"[INFO] Composition pie chart saved as '{save_path}'")
    plt.show()

# --------------------------------------------------
# 8. Density Report
# --------------------------------------------------
def print_density(structure):
    print(f"[INFO] Density = {structure.density:.3f} g/cm³")

# --------------------------------------------------
# Run All Analyses
# --------------------------------------------------
def run_analysis():
    struct = generate_lsmo_structure(x=0.3)
    export_to_lammps(struct)
    generate_lammps_input()

    try:
        run_lammps()
    except Exception as e:
        print(f"[ERROR] LAMMPS run failed: {e}")

    minimized = read_lammps_minimized()
    if minimized:
        print_density(minimized)
        plot_smoothed_xrd(minimized)
        plot_composition(minimized)
    else:
        print("[ERROR] Could not read minimized structure.")

# --------------------------------------------------
# Main
# --------------------------------------------------
if __name__ == "__main__":
    run_analysis()
