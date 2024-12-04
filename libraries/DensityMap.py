import numpy as np
import warnings
import random

import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter

from sklearn.cluster import KMeans

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import seaborn           as sns
sns.set_theme()

from sys import exit

# Suppress warning messages
warnings.filterwarnings("ignore", message="Failed to guess the mass for the following atom types")
warnings.filterwarnings("ignore", message="Guessed all Masses to 1.0")
warnings.filterwarnings("ignore", message="Reload offsets from trajectory\n ")
font_properties = FontProperties(family='serif', style='normal', weight='normal')

class DensityMap:
    def __init__(self, topofile = 'topo.gro', trajfile = 'trajectory.xtc', format = 'XDATCAR'
                 , timestep = 0.1, voxelsize = 0.2, atomname = 'Li', trajtimestep = 1, verbosity = 0, clustering = 0):
        """
        Initialize the density map and ion clustering analysis.

        Parameters:
        topo_file (str): topology file.
        traj_file (str): trajectory file.
        format (str): file format ('XDATCAR', 'XTC', 'TRR', 'LAMMPSDUMP').
        timestep (float): time step between frames.
        voxel_size (float): Size of the voxel for density calculation.
        atom_name (str or int): Name of the atom to be analyzed.
        trajtimestep (int): Frames interval to print the VoxelIndices.dat file.
        vervosity (int): Set verbosity level (0 or 1)
        clustering (str): Set whether clustering <di> of the trajectory is applied (0 = NO or 1 = YES)
        """
        self.topo_file = topofile
        self.traj_file = trajfile
        self.format = format 
        self.timestep = timestep
        self.voxel_size = voxelsize
        self.atom_name = atomname
        self.trajtimestep = trajtimestep        
        self.verbosity = verbosity
        self.clustering = clustering
        self.out_file = 'DensityMap.dat'

        self._density_matrix()
    
    def _read_XDATCAR(self):

        # Constants defining the information lines of the file
        scale_line  = 1  # Line for the scale of the simulation box
        s_cell_line = 2  # Start of the definition of the simulation box
        e_cell_line = 4  # End of the definition of the simulation box
        name_line   = 5  # Composition of the compound
        concentration_line = 6  # Concentration of the compound
        x_line = 7  # Start of the simulation      

        with open(self.format, 'r') as XDATCAR_file:
            XDATCAR_lines = XDATCAR_file.readlines()

        # Extracting the scale
        try:
            scale = float(XDATCAR_lines[scale_line])
        except ValueError:
            exit('Wrong definition of the scale in the XDATCAR.')

        # Extracting the cell
        try:
            box = np.array([line.split() for line in XDATCAR_lines[s_cell_line:e_cell_line+1]], dtype=float)
            box *= scale
        except ValueError:
            exit('Wrong definition of the cell in the XDATCAR.')

        # Extracting compounds and concentration
        compounds = XDATCAR_lines[name_line].split()
        concentration = np.array(XDATCAR_lines[concentration_line].split(), dtype=int)

        if len(compounds) != len(concentration):
            exit('Wrong definition of the composition of the compound in the XDATCAR.')

        n_ions = sum(concentration)

        # Extracting coordinates
        coordinates = np.array([line.split() for line in XDATCAR_lines[x_line:] if not line.split()[0][0].isalpha()], dtype=float)
        if not (len(coordinates) / n_ions).is_integer():
            exit('The number of lines is not correct in the XDATCAR file.')

        coordinates = coordinates.ravel().reshape((-1, n_ions, 3))

        # Creating a dictionary to map compounds to their concentrations and coordinates
        compound_data = {}
        start = 0
        for compound, num_atoms in zip(compounds, concentration):
            end = start + num_atoms
            compound_data[compound] = {
                'concentration': num_atoms,
                'coordinates': coordinates[:, start:end, :]
            }
            start = end

        # Returning coordinates for a specific element if requested
        if self.atom_name and self.atom_name in compound_data:
            coordinates =np.array(compound_data[self.atom_name]['coordinates'])
        
        for i in range(coordinates.shape[0]):
            coordinates[i]= np.dot(coordinates[i], box)
        
        return box, coordinates

    def _write_topo_file(self, coordinates, box):
        title="Topo File"
        with open(self.topo_file, "w") as file:
            file.write(f"{title}\n")
            total_atoms = len(coordinates)
            file.write(f"{total_atoms}\n")

            for atom_index in range(total_atoms):
                # Coordinates of the first frame
                x, y, z = coordinates[atom_index]/10
                file.write(f"{1:5d}{self.atom_name:>5}{self.atom_name:>5}{atom_index + 1:5d}"
                        f"{x:8.3f}{y:8.3f}{z:8.3f}\n")

            # Using cell dimensions from the 'cell' array for the box dimensions
            box_x, box_y, box_z = (np.diag(box))/10.0
            file.write(f"{box_x:10.5f}{box_y:10.5f}{box_z:10.5f}\n")

    def _write_xtc_file(self, coordinates, box):
  
        n_frames, n_atoms, _ = coordinates.shape
            
        # Extracting cell dimensions assuming an orthorhombic box
        box_dimensions = np.diag(box)  # Extract the diagonal to get the box dimensions
        
        # Creating a Universe with no topology, only using the coordinates
        universe = mda.Universe.empty(n_atoms, n_residues=n_atoms, atom_resindex=np.arange(n_atoms), trajectory=True)
        
        # Setting up the XTC writer
        with XTCWriter(self.traj_file, n_atoms) as xtc_writer:
            for frame in range(n_frames):
                # Updating the coordinates for each frame
                universe.atoms.positions = (coordinates[frame])
                # Use the same cell dimensions for all frames
                universe.trajectory.ts.dimensions = np.hstack([box_dimensions, [90.0, 90.0, 90.0]])  # Assuming orthorhombic box
                # Writing the frame
                xtc_writer.write(universe)
    
    def _density_map(self, u_traj, atom_id):        
        nx, ny, nz = (np.floor(u_traj.dimensions[:3] / self.voxel_size)).astype(int)
        dr = u_traj.dimensions[:3] / [nx, ny, nz]      
        
        m, n, k = nx, ny, nz
        n_frames = u_traj.trajectory.n_frames

        voxel_map = np.zeros((nx, ny, nz), dtype=int)
        box = u_traj.dimensions[0:3]

        index_voxmap_per_ion = []
        for frame_index, ts in enumerate(u_traj.trajectory):
            positions = u_traj.atoms.positions[atom_id]            

            for dim in range(3):
                over = positions[:, dim] >= box[dim] - dr[dim] / 2
                positions[over, dim] -= box[dim]

            voxel_indices = np.floor((positions + 0.5 * dr) / dr).astype(int)
            voxel_indices = np.clip(voxel_indices, 0, np.array([nx, ny, nz]) - 1)
            if frame_index % self.trajtimestep == 0:
                index_voxmap_per_ion.append(voxel_indices[:, 0] * n * k + voxel_indices[:, 1] * k + voxel_indices[:, 2] + 1)

            np.add.at(voxel_map, tuple(voxel_indices.T), 1)

        density_matrix = voxel_map/ (np.prod(dr) * n_frames)

        return nx, ny, nz, dr, n_frames, density_matrix, index_voxmap_per_ion


    def _write_output_density(self, density_matrix, out_dens_file, dr, targ = 0):

        m, n, k = density_matrix.shape

        with open(out_dens_file, 'w') as out:
            out.write(f'Density          x          y          z        index       voxel_ds: {dr[0]:.4f} {dr[1]:.4f} {dr[2]:.4f}\n')
            for i in range(m):
                for j in range(n):
                    for l in range(k):
                        index = i * n * k + j * k + l + 1  # Adjusted for 1-based indexing
                        value = density_matrix[i, j, l]
                        if targ == 0 or (targ == 1 and value > 0):
                            out.write(f'{value:^10.6e} {i + 1:^10} {j + 1:^10} {l + 1:^10}  {index:^10}\n')
  

    def _unwrap_traj(self, u_traj, atom_id, box):
        previous_positions = u_traj.atoms.positions[atom_id]
        coords = []; ids = []

        for ts in u_traj.trajectory:
            positions = u_traj.atoms.positions[atom_id]
            displacement = positions - previous_positions

            # To handle broadcasting, ensure both operands have the same shape
            wrapped_indices_pos = displacement > 0.5 * box[np.newaxis, :]
            wrapped_indices_neg = displacement < -0.5 * box[np.newaxis, :]

            correction_pos = wrapped_indices_pos * box[np.newaxis, :]
            correction_neg = wrapped_indices_neg * (-1) * box[np.newaxis, :]

            positions_unwrapped = positions - correction_pos - correction_neg
            coords.append(positions_unwrapped)
            ids.append(atom_id)
            previous_positions = positions_unwrapped.copy()

        coords_arr = np.array(coords)
        #ids = np.array(ids)
        mx, my, mz = coords_arr[:, :, 0], coords_arr[:, :, 1], coords_arr[:, :, 2]

        return ids, mx, my, mz
    
    def _feature_calc(self, mx, my, mz, n_atoms):
        tmax = len(mx)
        # Initialize the displacement array with zeros
        dr_avg = np.zeros(n_atoms)

        # Calculate the displacement for each time step and accumulate
        for i in range(1, tmax):
            ax = abs(np.subtract(mx[i][:],mx[0][:]))
            ay = abs(np.subtract(my[i][:],my[0][:]))
            az = abs(np.subtract(mz[i][:],mz[0][:]))
            ar = np.sqrt(ax**2 + ay**2 + az**2) / i
            dr_avg += ar
        return dr_avg
    
    def _clustering_di(self, features, ids):
        labels = []
        centroids = []
        kmeans = KMeans(n_clusters=3, n_init=50, random_state=123)
        for clust_crit in features:
            labels.append(kmeans.fit(clust_crit).predict(clust_crit))
            centroids.append(kmeans.fit(clust_crit).cluster_centers_)

        # it is important to check the labels of the min med and high 
        # finding the lables and id of atoms corresponding to the min med and high clusters
        arr_features = []; arr_features_sorted = []
        label_min = []; label_med = []; label_high =[]
        i = 0
        for feat in features:
            arr_features.append(np.array(list(zip(ids[0],feat[:,1],labels[i]))))    
            arr_features_aux = arr_features[i]
            arr_features_sorted.append(arr_features_aux[arr_features_aux[:,1].argsort()])    
            label_min.append(arr_features_sorted[i][0,2])
            label_high.append(arr_features_sorted[i][len(arr_features_sorted[i][:,0])-1, 2])
            known_values = [label_min[i], label_high[i]]
            label_med.append(3 - sum(known_values))
            i+=1

        #write files with the ids        
        output_info = ['di_avg.dat']        
        labels_sort = np.array(list(zip(label_min,label_med,label_high)))
        ids_per_cluster = []
        n = 0   
        for file_id in range(len(output_info)):            
            oinf = open(output_info[file_id],'w')                  
            for cluster_label in range(3):
                natoms_cluster = 0
                
                cluster_elements = []; cluster_idfl = []            
                for atom_n in range(len(arr_features[n][:,2])):
                    if arr_features[n][atom_n,2] == labels_sort[n,cluster_label]:
                        natoms_cluster +=1                    
                        cluster_elements.append(int(arr_features[n][atom_n,0]))
                        cluster_idfl.append(arr_features[n][atom_n])
                
                if file_id == 0:
                    ids_per_cluster.append(cluster_elements)
                

                oinf.write('{0} {1}'.format('#  id      di_avg     label       number of atoms = ', natoms_cluster))
                oinf.write("\n")
                for i in range(len(cluster_idfl)):
                    oinf.write('{:>6}  {:>.6f}  {:>1}'.format(int(cluster_idfl[i][0]),float(cluster_idfl[i][1]),int(cluster_idfl[i][2])))
                    oinf.write("\n")                    
            n+=1        
        oinf.close()

        ids_per_cluster_hm = ids_per_cluster[1] + ids_per_cluster[2]
        
        ids_per_cluster.append(ids_per_cluster_hm)
    
        return centroids, ids_per_cluster
    
    def _plot_clusters(self,centers):
        output_info = ['di_avg.dat']
        output_png = ['di_avg.png']
        plot_labelsy = ['$<d_i>$']
        colors = ['b', 'g', 'c','k', 'r', 'm', 'y']   
        scale_plot = 30
        
        for crit in range(len(output_info)):        
            fig, ax1 = plt.subplots(figsize=(8.0, 8.0))     
            with open(output_info[crit]) as f:
                dval_low = []; dval_med = []; dval_high = []                     
                l = 0; 
                for line in f:                
                    data = line.split()                                                         
                    if len(data)==9:      
                        l+=1                                   
                    if l == 1 and len(data)==3:
                        dval_low.append(float(data[1]))                                      
                    if l == 2 and len(data)==3:
                        dval_med.append(float(data[1]))                
                    if l == 3 and len(data)==3:
                        dval_high.append(float(data[1]))                   
                dval = [dval_low,dval_med,dval_high]
            
                c_code = [[0]*len(dval_low),[1]*len(dval_med),[2]*len(dval_high)] 

                total_dval = [i for subarray in dval for i in subarray]
                total_c = [i for subarray in c_code for i in subarray]
                total_elemnets = len(total_dval)

                x_center = [total_elemnets/2,total_elemnets/2,total_elemnets/2] 
            
                # Create a list of unique random integers
                unique_integers = random.sample(range(1, total_elemnets + 1), total_elemnets)
                for i in range(total_elemnets):
                    ax1.scatter(unique_integers[i],total_dval[i], color = colors[total_c[i]],s = scale_plot)
                            
                ax1.scatter(x_center[:],sorted(centers[crit][:,1]),marker = 'x',color = 'k',s = 200, linewidths = 4, zorder = 2.5)
    
            ax1.set_xlabel('{}'.format('Atom index'),fontproperties=font_properties, size=30)
            ax1.set_ylabel('{}'.format(plot_labelsy[crit]),fontproperties=font_properties, size=30)
            plt.xticks(fontproperties=font_properties, size=20)
            plt.yticks(fontproperties=font_properties, size=20)

            fig.tight_layout()
            plt.savefig('{}'.format(output_png[crit]))
            plt.close()


    def _write_chgcar(self, density_matrix, file_chgcar, box):        
        # Smoothing
        density_matrix = (density_matrix +
                    np.roll(density_matrix, shift=(1, 0, 0), axis=(0, 1, 2)) +
                    np.roll(density_matrix, shift=(-1, 0, 0), axis=(0, 1, 2)) +
                    np.roll(density_matrix, shift=(0, 1, 0), axis=(0, 1, 2)) +
                    np.roll(density_matrix, shift=(0, -1, 0), axis=(0, 1, 2)) +
                    np.roll(density_matrix, shift=(0, 0, 1), axis=(0, 1, 2)) +
                    np.roll(density_matrix, shift=(0, 0, -1), axis=(0, 1, 2)))/ 7
                    
        out_file = open(file_chgcar, 'w')
        count = 0
        out_file.write('Strucure\n 1.0\n {:.6} 0.0 0.0\n 0.0 {:.6} 0.0 \n 0.0 0.0 {:.6}\n    O\n1\n Direct\n 0.5 0.5 0.5\n\n'.
                        format(box[0],box[1],box[2]))
        out_file.write('{} {} {}\n'.format(int(density_matrix.shape[0]),int(density_matrix.shape[1]),int(density_matrix.shape[2])))

        for k in range(1,int(density_matrix.shape[2])+1):
            for j in range(1,int(density_matrix.shape[1])+1):
                for i in range(1,int(density_matrix.shape[0])+1):
                    count += 1
                    out_file.write('{:5.6e} '.format(density_matrix[i-1,j-1,k-1]))
                    if np.mod(count,5) == 0:
                        out_file.write('\n')
        out_file.close()

    def _write_output_voxel_visited(self, index_voxmap_per_ion, output_voxel_indices):

        # Convert to numpy array if not already and transpose
        index_voxmap_per_ion = np.array(index_voxmap_per_ion).T

        # Open the file for writing
        with open(output_voxel_indices, "w") as out_file:
            out_file.write('Indices of voxels visited by each atom over time\n')
            # Write each row as a tab-separated string
            for row in index_voxmap_per_ion:
                out_file.write("\t".join(map(str, row)) + "\n")


    def _write_output_cluster(self, file_chgcar, voxels_visited_file, ids_per_cluster, u_traj, box):
        print('>> Calculating density map for <di> clusters <<')        

        # Define file names based on verbosity
        if self.verbosity == 1:
            cluster_labels = ['L', 'M', 'H', 'HM']
        else:
            cluster_labels = ['L', 'HM']

        density_files_di_avg = [f'{label}-{self.out_file}' for label in cluster_labels]
        voxel_index_di_avg = [f'{label}-{voxels_visited_file}' for label in cluster_labels]
        chgcar_files_di_avg = [f'{file_chgcar}-{label}' for label in cluster_labels]

        # Update ids_di_avg to pick the right indices based on cluster_labels
        all_cluster_labels = ['L', 'M', 'H', 'HM']
        ids_di_avg = [ids_per_cluster[all_cluster_labels.index(label)] for label in cluster_labels]
        print(f'Time interval between frames DensityMaps: {self.timestep:.4f}')
        print(f'Time interval between frames VoxelIndices: {self.trajtimestep*self.timestep:.4f}')
        print(f'Number of frames VoxelIndices: {self.n_frames // self.trajtimestep:d}')
        # Process and write files for each cluster
        for i, cluster_id in enumerate(ids_di_avg):
            ids_di_avg_corrected = np.array(cluster_id) - 1,            
            print(f'Writing files: {density_files_di_avg[i]}, {voxel_index_di_avg[i]}, {chgcar_files_di_avg[i]}' )

            _, _, _, dr, _, density_matrix, index_voxmap_per_ion = self._density_map(u_traj, ids_di_avg_corrected)
            
            self._write_output_voxel_visited(index_voxmap_per_ion, voxel_index_di_avg[i])
            self._write_output_density(density_matrix, density_files_di_avg[i], dr, targ = 1)            
            self._write_chgcar(density_matrix, chgcar_files_di_avg[i], box)

    def _density_matrix(self):
        
        print('>> Reading the trajectory <<')
        supported_formats = ['XTC', 'TRR', 'LAMMPSDUMP', 'XDATCAR']
        if self.format not in supported_formats:
            exit(f"Error: Unsupported trajectory format '{self.format}'.")
            

        try:
            if self.format in ['XTC', 'TRR']:
                if isinstance(self.atom_name, str) and self.atom_name.isalpha(): #self.atom_name.isalpha():
                    atom_type = 'name'
                else:
                    atom_type = 'type'            
                u_traj = mda.Universe(self.topo_file, self.traj_file, format=f'{self.format}')
                atom_id = u_traj.select_atoms(f'{atom_type} {self.atom_name}').indices                
                total_atoms = len(atom_id)
            elif self.format == 'LAMMPSDUMP':
                if self.atom_name.isalpha():
                    atom_type = 'name'
                else:
                    atom_type = 'type'
                u_traj = mda.Universe(self.topo_file, self.traj_file, format=f"{self.format}")
                atom_id = u_traj.select_atoms(f'{atom_type} {self.atom_name}').indices
                total_atoms = len(atom_id)                             
            elif self.format == 'XDATCAR':                
                box, coordinates = self._read_XDATCAR()                
                self._write_topo_file(coordinates[0], box)
                self._write_xtc_file(coordinates, box)            
                u_traj = mda.Universe(self.topo_file, self.traj_file, format="XTC")
                atom_id = u_traj.select_atoms(f'name {self.atom_name}').indices
                total_atoms = len(atom_id)
        except Exception as e:
            exit(f"Error loading trajectory data: {e}")
             
        file_chgcar = 'CHGCAR_VESTA'
        voxel_index_file  = 'VoxelIndices.dat'
        
        total_atoms = len(atom_id)

        print('>> Calculating density map <<')
        box = u_traj.dimensions[:3]        
        nx, ny, nz, dr, self.n_frames, density_matrix, index_voxmap_per_ion = self._density_map(u_traj,atom_id)
        
        if self.clustering == 0:            
            print(f'Format: {self.format}')
            print(f'Number frames DensityMap: {self.n_frames}')
            print(f'Diffusive atoms: {total_atoms}')
            print(f'Voxel size: {self.voxel_size}')
            print(f'Number voxels: {(nx)*(ny)*(nz)}\n')

            out_dens_file = self.out_file
            print(f'Writing file: {self.out_file}')
            print(f'Time interval between frames DensityMap: {self.timestep:.4f}')
            self._write_output_density(density_matrix, out_dens_file, dr, targ = 1)           
            print(f'\nWriting file: {voxel_index_file}')
            print(f'Time interval between frames VoxelIndices: {self.trajtimestep * self.timestep:.4f}')
            print(f'Number of frames VoxelIndices: {self.n_frames // self.trajtimestep:d}')            
            self._write_output_voxel_visited(index_voxmap_per_ion, voxel_index_file) 
            print(f'\nWriting file: {file_chgcar}')
            self._write_chgcar(density_matrix, file_chgcar, box)
            print('>> Done <<\n')
        else:
            print(f'Format: {self.format}')
            print(f'Number frames DensityMap: {self.n_frames}')
            print(f'Diffusive atoms: {total_atoms}')
            print(f'Voxel size: {self.voxel_size}')
            print(f'Number voxels: {(nx)*(ny)*(nz)}\n')

            print('>> Computing <di> clustering <<')
            ids, mx, my, mz = self._unwrap_traj(u_traj, atom_id, box)
            frame_i, frame_f = int(self.n_frames * 0.1), int(self.n_frames * 0.9)        
            print(f'Initial frame: {frame_i}')
            print(f'Final frame: {frame_f}')         
            print(f'Total Steps used: {len(mx[frame_i:frame_f])}')
            print('>> Done <<\n')
            
            # Calculate average displacements <di>        
            ids = np.array(ids) + 1  # Increase Id index +1 because MDanalysis starts indices at 0
            di_avg = self._feature_calc(mx[frame_i:frame_f,:], my[frame_i:frame_f,:], mz[frame_i:frame_f,:], total_atoms)

            Xdi_avg = np.array(list(zip(di_avg,di_avg)))
            features = [Xdi_avg]
            centers, ids_per_cluster = self._clustering_di(features, ids)
            self._plot_clusters(centers)
            
            cluster_labels = ['L', 'M', 'H', 'HM']

            hm_count = len(ids_per_cluster[1]) + len(ids_per_cluster[2])

            output_lines = [f'Number of atoms per <di> cluster:']
            
            for i, label in enumerate(cluster_labels[:-1]):
                output_lines.append(f'{label}: {len(ids_per_cluster[i])}')
            output_lines.append(f'{cluster_labels[-1]}: {hm_count}')
            # Print the formatted output
            print('\n'.join(output_lines))
                    
            # Write output files for clusters
            self._write_output_cluster(file_chgcar, voxel_index_file, ids_per_cluster, u_traj, box)

            print('>> Done <<\n')
