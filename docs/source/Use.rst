Getting started
===================

After performing an atomistic simulation, **CrySF** can be used for analysis via the following steps:

1. **Prepare the Simulation Data:**

   a. Store the atomistic trajectory in a dedicated folder. Include the following files based on the simulation software used:

      - `XDATCAR` from VASP.
      - For LAMMPS or GROMACS, include the topology and trajectory files.

2. **Set Up the Analysis Scripts:**

   - Copy the files `densitymap.py` and `crysf.py` into the same folder as your simulation data.

3. **Execute `DensityMap` to Create a Density Map:**

   - Run the following command to execute `densitymap.py`:

     ```
     python densitymap.py -to 'str' -tr 'str' -f 'str' -ts float -v float -a int or 'str' -tts int -verb int -clus int
     ```

4. **Execute `CrySF` to Analyze the Density Map:**

   - After generating the density map, run `crysf.py` with the command:

     ```
     python crysf.py -nts float -minv float -maxv float -clus int -dop int -verb int -deltf int -scaler 'str'
     ```

The settings to use the `DensityMap` and `CrySF` commands are detailed in the **Settings** section.

The time it takes to analyze a trajectory depends on the simulation's number of atoms and the frequency with which the trajectory is printed. Usually, the analysis takes a few minutes.
