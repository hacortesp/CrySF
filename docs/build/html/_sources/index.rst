CrySF
=================================
Ionic diffusion is fundamental to materials science, driving innovations in energy storage, catalysis, and solid-state devices. Understanding atomic motion is essential for optimizing properties like ionic conductivity, stability, and performance, especially in technologies such as solid-state batteries and fuel cells.

We developed `CrySF` (Crystallography Sites Finder), a Python-based tool designed to analyze molecular dynamics trajectories of crystalline ionic conductors. CrySF identifies crystallographic sites occupied by diffusive atoms and tracks both individual and collective transitions, offering critical insights into how structures influence atoms individual and cooperative diffusion.

Using a two-step methodology— `DensityMap` and `CrySF` —the software automatically generates a density map and analyzes it to extract the following key details:

- The total number of crystallographic sites visited by diffusive atoms.
- Classification of sites based on their geometry and volume.
- Calculation of site occupancies.
- Identification of coordination relationships between sites.
- The number of distinct sites visited by each diffusive atom.
- Quantification of forward and backward jumps for each atom.
- Detection of simultaneous jumps by multiple atoms.
- Analysis of cooperative jump strings, including their frequency and lengths.

By combining structural characterization with atomic transport dynamics, CrySF provides an automated and detailed framework to accelerate the design of materials with enhanced ionic transport properties, advancing innovation in energy storage and solid-state technologies.


Requirements
------------

- Python (>=3.8): `<https://www.python.org/downloads/>`_
- MDAnalysis: `<https://www.mdanalysis.org/>`_
- scikit-learn: `<https://scikit-learn.org/stable/>`_
- Matplotlib: `<https://matplotlib.org/>`_
- SciPy: `<https://scipy.org/>`_
- Seaborn: `<https://seaborn.pydata.org/>`_

Installation
------------

CrySF can be installed from source::

    git clone https://gitlab.bcamath.org/hacortes/crysf
    cd crysf
    pip install .

or used directly from source without explicit installation::

    git clone https://gitlab.bcamath.org/hacortes/crysf
    cd crysf
    pip install -r docs/requirements.txt


Execution
---------

An *ab initio* MD (VASP)  and clasic MD (LAMMPS) simulations of the cubic phase of Li\ :sub:`7`\ La\ :sub:`3`\ Zr\ :sub:`2`\ O\ :sub:`12`\  (cLLZO) at 300K and Li\ :sub:`10`\ GeP\ :sub:`2`\ S\ :sub:`12`\  (LGPS) at 650K  are provided as example simulations. More details of the examples can be found in `examples <https://gitlab.bcamath.org/hacortes/crysf/-/tree/main/examples>`_.

Input trajectories
------------------

`CrySF` soports different trajectory formats: XTC, TRR, LAMMPSDUMP and XDATCAR. 

Authors
-------

CrySF is being developed by:

- Dr. Henry Andres Cortes, BCAM - Basque Center for Applied Mathematics, Spain
- Dr. Mauricio Rincon Bonilla, BCAM - Basque Center for Applied Mathematics, Spain
- Prof. Elena Akhmatskaya, BCAM - Basque Center for Applied Mathematics, Spain


Contact
-----------------------------------

If you have questions, please don't hesitate to reach out at: hacortes@bcamath.org GitLab: `CrySF GitLab Repository <https://gitlab.bcamath.org/hacortes/crysf>`_ and Github: `CrySF Github Repository <https://github.com/hacortesp/CrySF>`_

.. toctree::
   :hidden:
   :caption: Usage
   :maxdepth: 4

   source/Use
   source/Settings
   source/CrySF
