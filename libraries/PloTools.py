#Work in progress
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MaxNLocator, LogLocator, FuncFormatter

import seaborn as sns
sns.set_theme()

rc("text", usetex=False)
rc("font", family="serif")

font_properties = FontProperties(family='serif', style='normal', weight='normal')
LABEL_FONTSIZE = 30
TICK_FONTSIZE = 20

def log_format(y, _):
    return f'$10^{{{int(np.log10(y))}}}$' if y != 0 else '0'

def plot_cluster_shapes(variances, labels_shapes, output_filename='site_shapes_variance.png'):
    """
    Plot scatter plots comparing shape variances of clusters across axes.

    Parameters:
        variances (np.ndarray): 2D array where each row represents [lambda_x, lambda_y, lambda_z] for a cluster.
        labels_shapes (array-like): Labels used to color each cluster based on its shape group.
        output_filename (str): Filename for saving the resulting figure.

    Output:
        Saves a two-panel scatter plot to a PNG file.
    """
    title = r'Site shapes variance'
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.0))

    scatter1 = axes[0].scatter(
        variances[:, 0],
        variances[:, 2],
        c=labels_shapes,
        cmap='rainbow',
        edgecolor='k',
        s=40
    )
    axes[0].set_title(r'$\lambda_x$ vs $\lambda_z$', fontproperties=font_properties, fontsize=LABEL_FONTSIZE)
    axes[0].set_xlabel(r'$\lambda_x$', fontproperties=font_properties, fontsize=LABEL_FONTSIZE, labelpad=10)
    axes[0].set_ylabel(r'$\lambda_z$', fontproperties=font_properties, fontsize=LABEL_FONTSIZE, labelpad=10)
    axes[0].tick_params(axis='both', labelsize=TICK_FONTSIZE, labelrotation=0)

    scatter2 = axes[1].scatter(
        variances[:, 1],
        variances[:, 2],
        c=labels_shapes,
        cmap='rainbow',
        edgecolor='k',
        s=40
    )
    axes[1].set_title(r'$\lambda_y$ vs $\lambda_z$', fontproperties=font_properties, fontsize=LABEL_FONTSIZE)
    axes[1].set_xlabel(r'$\lambda_y$', fontproperties=font_properties, fontsize=LABEL_FONTSIZE, labelpad=10)
    axes[1].set_ylabel(r'$\lambda_z$', fontproperties=font_properties, fontsize=LABEL_FONTSIZE, labelpad=10)
    axes[1].tick_params(axis='both', labelsize=TICK_FONTSIZE, labelrotation=0)

    fig.suptitle(title, fontproperties=font_properties, fontsize=LABEL_FONTSIZE + 10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    plt.savefig(output_filename)
    plt.close()


def plot_sites_visited(sites_visited, output_filename='sites_visited.png'):
    """
    Plot a histogram showing how many atoms visited how many sites.

    Parameters:
        sites_visited (array-like): Each entry indicates how many sites a single atom has visited.
        output_filename (str): Filename for saving the resulting bar plot.

    Output:
        Saves a histogram bar plot to a PNG file.
    """
    unique_values, counts = np.unique(sites_visited, return_counts=True)
    fig, ax = plt.subplots(figsize=(8.0, 8.0))

    ax.bar(unique_values, counts, color='b', alpha=0.7, edgecolor='k', width=0.8)
    ax.set_xlabel('Number of sites visited', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Number of Atoms', fontproperties=font_properties, fontsize=LABEL_FONTSIZE)
    ax.set_xlim([min(unique_values) - 1, max(unique_values) + 1])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE, labelrotation=0)

    fig.tight_layout()
    plt.savefig(output_filename)
    plt.close()


def plot_cluster_amplitudes(amplitud_values, labels_shapes, output_filename='site_amplitudes.png'):
    """
    Plot site amplitudes as a scatter plot with clusters color-coded.

    Parameters:
        amplitud_values (array-like): Amplitude values for each site.
        labels_shapes (array-like): Labels representing the shape cluster each site belongs to.
        output_filename (str): Filename for saving the resulting scatter plot.

    Output:
        Saves a scatter plot to a PNG file.
    """
    amplitud_indices = range(len(amplitud_values))
    fig, ax = plt.subplots(figsize=(8.0, 8.0))

    ax.scatter(
        amplitud_indices,
        amplitud_values,
        c=labels_shapes,
        cmap='rainbow', edgecolor='k', s=40
    )

    ax.set_xlabel('Site Index', fontproperties=font_properties, fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Site amplitude (A)', fontproperties=font_properties, fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE, labelrotation=0)

    fig.tight_layout()
    plt.savefig(output_filename)
    plt.close()


def plot_simultaneous_jumps(filename: str, delt_ps: float):
    """
    Plot histogram of the number of simultaneous jumps over time.

    Parameters:
        filename (str): Path to the data file containing frame and NSjumps.
        delt_ps (float): Time step size in picoseconds per frame.

    Output:
        Saves a histogram with error bars to 'simultaneous_jumps_histogram.png'.
    """
    data = pd.read_csv(filename, sep=r'\s+', comment='#', header=None, names=['Frame', 'NSjumps'])
    data['Time_ps'] = data['Frame'] * delt_ps

    total_time = data['Time_ps'].max()
    num_bins = 10
    bins = np.linspace(0, total_time, num_bins + 1)
    bin_width = bins[1] - bins[0]
  
    print(f'bins width in simultaneous jumps: {bin_width:.4f} ps')
    data['Bin'] = pd.cut(data['Time_ps'], bins)
    grouped = data.groupby('Bin', observed=False)['NSjumps']

    bin_centers = [interval.left + bin_width / 2 for interval in grouped.groups.keys()]
    means = grouped.mean().values
    stds = grouped.std().fillna(0).values

    plt.figure(figsize=(8.0, 8.0))
    plt.bar(bin_centers, means, width=bin_width, color='b', alpha=0.7, edgecolor='black', yerr=stds, capsize=3)
    plt.xlabel(r'$t$ [ps]', fontsize=LABEL_FONTSIZE)
    plt.ylabel(r'Number of jupms', fontsize=LABEL_FONTSIZE)
    plt.tick_params(axis='both', labelsize=TICK_FONTSIZE)
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig('simultaneous_jumps.png')
    plt.close()

def plot_string_frequency(filename: str):
    """
    Plot the frequency distribution of string lengths on a log scale.

    Parameters:
        filename (str): Path to the file containing string length and probability data.

    Output:
        Saves a log-scaled probability plot to 'string_prob.png'.
    """
    data = pd.read_csv(filename, sep=r'\s+', comment='#', header=None,
                       names=['String Length', 'f_n(ps^-1)'])
    plt.figure(figsize=(10.0, 8.0))
    plt.scatter(data['String Length'], data['f_n(ps^-1)'], color='black', s=150)
    plt.plot(data['String Length'], data['f_n(ps^-1)'], marker='o', linestyle='-', color='black')

    max_length = data['String Length'].max()
    plt.xlim(0, max_length + 1)
    plt.xticks(np.arange(0, max_length + 2, 1))

    ymin = max(data['f_n(ps^-1)'].min() * 0.8, 1e-6)
    ymax = data['f_n(ps^-1)'].max() * 1.2
    plt.ylim(ymin, ymax)
    plt.yscale('log')

    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter())
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
    plt.gca().yaxis.set_major_formatter(FuncFormatter(log_format))

    plt.xlabel(r'n', fontproperties=font_properties, fontsize=LABEL_FONTSIZE)
    plt.ylabel(r'$f_{n}$ (ps$^{-1}$) ', fontsize=LABEL_FONTSIZE, fontproperties=font_properties)
    plt.tick_params(axis='both', labelsize=TICK_FONTSIZE)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    plt.savefig('string_freq.png')
    plt.close()

def plot_jump_coefficient(filename):
    """
    Plot self-correlation values per atom, color-coded by a flag (-1 means red).

    Parameters:
        filename (str): Path to the file with columns Atom, NJumps, RNJumps, SelfCorr.

    Output:
        Saves a scatter plot to 'self_corr.png'.
    """
    data = pd.read_csv(filename, sep=r'\s+', comment='#', header=None,
                       names=['Atom', 'NJumps', 'RNJumps', 'SelfCorr'])

    plt.figure(figsize=(8.0, 8.0))
    colors = ['red' if x == -1 else 'black' for x in data['SelfCorr']]
    edgecolors = ['red' if x == -1 else 'white' for x in data['SelfCorr']]

    for color, edgecolor, (x, y) in zip(colors, edgecolors, data[['Atom', 'SelfCorr']].values):
        plt.scatter(x, y, color=color, edgecolor=edgecolor, marker='o', s=150)

    plt.ylim(0, 1)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.xlim(0, data['Atom'].max() + 1)

    plt.xlabel(r'Atoms', fontproperties=font_properties, fontsize=LABEL_FONTSIZE)
    plt.ylabel(r'$\alpha_{k}$', fontproperties=font_properties, fontsize=LABEL_FONTSIZE)
    plt.tick_params(axis='both', labelsize=TICK_FONTSIZE)
    plt.savefig('jump_coefficient.png')
    plt.close()


#-----> Write files 

def write_occupancy_ovito(box_lengths):
    """
    Parse SitesMap.dat and write a LAMMPS trajectory file (occ.lammpstrj).

    Parameters:
        box_lengths (list or tuple of float): Lengths of the simulation box in x, y, z directions.

    Input:
        'SitesMap.dat' - Input file containing atom positions and occupancy data.

    Output:
        'occ.lammpstrj' - LAMMPS trajectory file including occupancy field.
    """
    if len(box_lengths) != 3:
        raise ValueError("box_lengths must be a list or tuple of three floats: [x_len, y_len, z_len]")

    box_bounds = [(0.0, box_lengths[0]), (0.0, box_lengths[1]), (0.0, box_lengths[2])]

    sitesmap_file = 'SitesMap.dat'
    output_file = 'occ.lammpstrj'

    atoms = []
    with open(sitesmap_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.split()
                if not parts[0].isdigit():
                    continue
                atom_id = int(parts[0])
                atom_type = int(parts[4])
                x, y, z = map(float, parts[1:4])
                occupancy = float(parts[6])
                atoms.append((atom_id, atom_type, x, y, z, 1.0, occupancy))

    with open(output_file, 'w') as f:
        f.write("ITEM: TIMESTEP\n")
        f.write("0\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{len(atoms)}\n")
        f.write("ITEM: BOX BOUNDS pp pp ss\n")
        for lo, hi in box_bounds:
            f.write(f"{lo:.16e} {hi:.16e}\n")
        f.write("ITEM: ATOMS id type x y z q occ\n")
        for atom in atoms:
            f.write(f"{atom[0]} {atom[1]} {atom[2]:.5f} {atom[3]:.5f} {atom[4]:.5f} {atom[5]:.6f} {atom[6]}\n")

    print(f"Written: {output_file}")



def diffusion_pathway_ovito(box_lengths):
    """
    Generate a LAMMPS data file describing diffusion pathways.

    Parameters:
        box_lengths (tuple of float): Lengths of the simulation box in x, y, z directions.

    Input:
        - 'SitesMap.dat': Atom site information.
        - 'jumping_path.dat': Bonding path (jump) information.

    Output:
        - 'jumping_path_ovito.dat': LAMMPS data file for OVITO visualization.
    """
    sites_file='SitesMap.dat'
    paths_file='jumping_path.dat'
    out_file='jumping_path_ovito.dat'
    CHARGE = 0.0
    MOL_ID = 1

    def parse_sites(file_path):
        atoms = []
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split()
                atom_id = int(parts[0])
                x, y, z = map(float, parts[1:4])
                atom_type = int(parts[4]) + 1
                atoms.append((atom_id, MOL_ID, atom_type, CHARGE, x, y, z))
        return atoms

    def parse_bonds(file_path):
        bonds = []
        bond_types = []
        bond_id = 1
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                i, j, occ = map(int, line.split())
                bonds.append((bond_id, occ, i, j))
                bond_types.append(occ)
                bond_id += 1
        return bonds, bond_types

    def write_lammps_data(outfile, atoms, bonds, bond_types, box_lengths):
        xlen, ylen, zlen = box_lengths
        fmin = 0.0
        with open(outfile, 'w') as f:
            f.write("# LAMMPS data file generated by script\n")
            f.write(f"{len(atoms)} atoms\n")
            f.write(f"{len(bonds)} bonds\n")
            f.write(f"{max(a[2] for a in atoms)} atom types\n")
            f.write(f"{max(bond_types)} bond types\n")
            f.write(f"{fmin} {xlen} xlo xhi\n")
            f.write(f"{fmin} {ylen} ylo yhi\n")
            f.write(f"{fmin} {zlen} zlo zhi\n\n")

            f.write("Masses\n\n")
            for t in range(1, max(a[2] for a in atoms)+1):
                f.write(f"{t} 0.0  # Type{t}\n")

            f.write("\nAtoms  # full\n\n")
            for a in atoms:
                f.write(f"{a[0]} {a[1]} {a[2]} {a[3]} {a[4]} {a[5]} {a[6]}\n")

            f.write("\nBonds\n\n")
            for b in bonds:
                f.write(f"{b[0]} {b[1]} {b[2]} {b[3]}\n")

    atoms = parse_sites(sites_file)
    bonds, bond_types = parse_bonds(paths_file)
    write_lammps_data(out_file, atoms, bonds, bond_types, box_lengths)
    print(f"Written: {out_file}")
