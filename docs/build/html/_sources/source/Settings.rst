Settings
========

Command-line Settings for `DensityMap`
--------------------------------------

The following table provides an overview of the command-line parameters available for the `DensityMap` script:

+-------------+---------------------------------------------------------------+-------------------+
| Parameter   | Description                                                   | Default           |
+=============+===============================================================+===================+
| `-to`       | Topology file name (.gro or .data formats). XDATCAR does not  | `topo.gro`        |
|             | require to include `-to`.                                     |                   |
+-------------+---------------------------------------------------------------+-------------------+
| `-tr`       | Trajectory file name.                                         | `trajectory.xtc`  |
+-------------+---------------------------------------------------------------+-------------------+
| `-f`        | Format options: XTC, TRR, LAMMPSDUMP, and XDATCAR.            | `XDATCAR`         |
+-------------+---------------------------------------------------------------+-------------------+
| `-ts`       | Time interval between frames in the trajectory.               | `0.1`             |
+-------------+---------------------------------------------------------------+-------------------+
| `-v`        | Voxel side size.                                              | `0.2`             |
+-------------+---------------------------------------------------------------+-------------------+
| `-a`        | Atom selection: integer (index) or string (name) based on     | `Li`              |
|             | trajectory format.                                            |                   |
+-------------+---------------------------------------------------------------+-------------------+
| `-tts`      | Frames interval to print the `VoxelIndices.dat`.              | `1`               |
+-------------+---------------------------------------------------------------+-------------------+
| `-clus`     | Set whether clustering <d₁> of the trajectory is applied      | `0`               |
|             | (0 = NO or 1 = YES).                                          |                   |
+-------------+---------------------------------------------------------------+-------------------+
| `-verb`     | Set verbosity level (0 = L and HM or 1 = L, M, H, and HM).    | `0`               |
+-------------+---------------------------------------------------------------+-------------------+

**Notes:**

1. The frequency at which the position of each diffusive atom is recorded in terms of voxel indices can be adjusted. A value equivalent to 1 ps is suggested (i.e., `-ts` * `-tts` ≈ 1 ps) since it is in the range of a typical phonon frequency [Lim-2018, Kam-2023].
2. L, M, and H represent the <d₁> clusters low, medium, and high, respectively [`Cortes-2024 <https://pubs.rsc.org/en/content/articlelanding/2024/ta/d3ta07036k>`_].



Command-line Settings for `CrySF`
---------------------------------

The following table provides an overview of the command-line parameters available for the `CrySF` script:

+-------------+---------------------------------------------------------------+-------------------+
| Parameter   | Description                                                   | Default           |
+=============+===============================================================+===================+
| `-nts`      | New time interval between frames in `VoxelIndices.dat`.       | `0.1`             |
+-------------+---------------------------------------------------------------+-------------------+
| `-minv`     | Minimum volume condition for the site.                        | `0.28`            |
+-------------+---------------------------------------------------------------+-------------------+
| `-maxv`     | Maximum volume condition for the site.                        | `3.05`            |
+-------------+---------------------------------------------------------------+-------------------+
| `-clus`     | Set whether clustering <d₁> of the trajectory is used         | `0`               |
|             | (0 = NO or 1 = YES).                                          |                   |
+-------------+---------------------------------------------------------------+-------------------+
| `-dop`      | Set whether to use a doped density map (0 = NO or 1 = YES).   | `0`               |
+-------------+---------------------------------------------------------------+-------------------+
| `-verb`     | Set verbosity level, 0 = L and HM or 1 = L, M, H, and HM.     | `0`               |
+-------------+---------------------------------------------------------------+-------------------+
| `-deltf`    | Frames intervals to calculate the simultaneous jumps.         | `1`               |
+-------------+---------------------------------------------------------------+-------------------+
| `-scaler`   | Scaler method used for site type identification:              | `MinMaxScaler`    |
|             | MinMaxScaler or StandardScaler.                               |                   |
+-------------+---------------------------------------------------------------+-------------------+

**Notes:**

1. The `-nts` value is an output of the `DensityMap` routine in the Voxel indices information (see `Time interval between frames`).
2. Further details about clustering are discussed in the **<d₁> Clustering** section.
3. Additional information about doped systems can be found in the **Doped Systems** section.
4. Given that simultaneous jumps are rare events, it is often advantageous to segment the trajectory into chunks for analysis. This setting is an alternative to `-tts > 1` in `DensityMap`. Using `-tts = 1` (all frames), you can split the trajectory into chunks and establish whether a jump occurred in each.

---

Units of Measurement
--------------------

The software uses the following units of measurement:

- **Distance**: Angstroms (`Å`)
- **Time**: Picoseconds (`ps`)
