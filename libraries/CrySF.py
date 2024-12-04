import numpy as np
import sys

from scipy.spatial import cKDTree

from collections import deque, Counter,defaultdict

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors 
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import MaxNLocator
import seaborn           as sns

sns.set_theme()

kmeans_kwargs   = dict(init='random', n_init=10, max_iter=300, tol=1e-04, random_state=0)

class CrySF:
    def __init__(self, timestep = 0.2, minv = 0.28, maxv = 3.05, 
                 verbosity = 0, deltframes = 1, scaler = 'MinMaxScaler',
                 clustering = 0,doped = 0):
     
        self.density_map_targ = 'DensityMap.dat'
        self.out_file = 'SitesMap.dat'            
        self.time_frames = timestep
        self.minv = minv
        self.maxv = maxv
        self.verbosity = verbosity
        self.clustering = clustering
        self.doped = doped
        self.deltf = deltframes
        self.scaler = scaler
        self.font_properties = FontProperties(family='serif', style='normal', weight='normal')

        self._site_analysis()

    def _read_density_map(self, file_density):        
        with open(file_density) as f:            
            first_line = next(f).split()
            dl = [float(first_line[6]), float(first_line[7]), float(first_line[8])]

            density, xyz, index_voxels = zip(*[(float(data[0]), [int(s) for s in data[1:4]], int(data[4]))
                                        for data in (line.split() for line in f)])
        self.dl = np.array(dl)
        self.voxel_vol = self.dl[0]*self.dl[1]*self.dl[2]
        xyz_map = np.array(list(map(list, xyz)))
        self.pbc = np.array([max(xyz_map[:,i]) for i in range(3)])        
        self.eps = np.sqrt(3)
        self.total_time = (self.n_frames-1) * (self.time_frames)
        density_values = np.array(list(density))        
        index_voxels = np.array(index_voxels)

        return density_values, xyz_map, index_voxels
    
    def _clustering_pbc(self, X):
        min_samples = 12 #the number of NN voxels 
        X_mod = np.mod(X, self.pbc)
        tree = cKDTree(data=X_mod, boxsize=self.pbc)
        indices = tree.query_ball_tree(tree, r = self.eps)
        
        labels = np.full(shape=X.shape[0], fill_value=-1, dtype=int)  # Initialize all points as noise
        visited = set()  # Keep track of visited points to avoid redundancy
        cluster_id = 0

        for i, neighbors in enumerate(indices):
            if len(neighbors) >= min_samples and i not in visited:  # Check if not visited for efficiency
                labels[i] = cluster_id
                visited.add(i)
                points_to_visit = deque(neighbors)  # Use deque for efficient queue operations
                
                while points_to_visit:
                    point = points_to_visit.popleft()
                    if point not in visited:  # Check if not visited
                        visited.add(point)
                        labels[point] = cluster_id
                        new_neighbors = tree.query_ball_point(X_mod[point], r = self.eps)
                        if len(new_neighbors) >= min_samples:
                            points_to_visit.extend(new_neighbors)
                cluster_id += 1
        return labels
    
    def _sites_cova_pca(self, voxel_coor, labels_cut):
        max_label = np.max(labels_cut)
        cluster_spatial_cova = []

        half_pbc = 0.5 * self.pbc
        pbc_correction_factors = np.array([half_pbc, self.pbc, -self.pbc])

        for label in range(max_label + 1):            
            cluster_coor = voxel_coor[labels_cut == label]            
            if len(cluster_coor) > 1:
               
                distances = cluster_coor[:, np.newaxis, :] - cluster_coor[np.newaxis, :, :]
                distance_pos = (distances > half_pbc) * pbc_correction_factors[1]
                distance_neg = (distances < -half_pbc) * pbc_correction_factors[2]
                adjusted_distances = distances - (distance_pos + distance_neg)

                cluster_coor_adjusted = cluster_coor[0] + adjusted_distances[:, 0, :] * self.dl
                # Check if there are enough samples and features
                n_samples, n_features = cluster_coor_adjusted.shape
                if n_samples >= 3 and n_features >= 3:
                    pca = PCA(n_components=3, svd_solver='full')
                    pca.fit(cluster_coor_adjusted)
                    cluster_spatial_cova.append(pca.explained_variance_[:3].tolist())
                else:
                    # Handle cases where there are not enough samples or features
                    # For example, append zeros, handle differently, or log an error                   
                    cluster_spatial_cova.append([0.0, 0.0, 0.0])  # Example fallback
      
        return np.array(cluster_spatial_cova)
 
    def _write_cluster_map(self, coordinates, labels, out_file):
        #out_file = 'SitesMap_OVITO.dat'
        total_atoms = len(coordinates)
        atom_types = len(set(labels))
        bounds_template = "{0:.8f} {1} {2}\n"
        atom_line_template = "{:<10d}{:<10d}{:<12.6f}{:<12.6f}{:<12.6f}\n"
        

        with open(out_file, 'w') as file:
            file.write("LAMMPS data file via OVITO\n\n")
            file.write(f"{total_atoms} atoms\n")
            file.write(f"{atom_types} atom types\n\n")
            file.write(bounds_template.format(0.0, self.pbc[0]*self.dl[0], "xlo xhi"))
            file.write(bounds_template.format(0.0, self.pbc[1]*self.dl[0], "ylo yhi"))
            file.write(bounds_template.format(0.0, self.pbc[2]*self.dl[0], "zlo zhi"))
            file.write("\nMasses\n\n")

            for i in range(atom_types):
                file.write(f"{i+1} 1.0\n")

            file.write("\nAtoms # atomic\n\n")

            for i, coord in enumerate(coordinates):
                scaled_coord = coord * self.dl[0]              
                file.write(atom_line_template.format(i + 1, labels[i] + 1, *scaled_coord))

    def _find_density_cutoff(self):
        self.density_values, self.xyz_map, self.index_voxels = self._read_density_map(self.density_map_targ)
        print('>> Finding the density cutoff for clustering <<')
        n, bins, _ = plt.hist(self.density_values, bins=30, alpha=0.7, color='blue')        
        max_bin_index = np.argmax(n)
        value_i = bins[max_bin_index]
        value_s = bins[max_bin_index + 8]
        
        self.cutoff_values = np.arange(value_i, value_s + value_i, value_i) 
        plt.close()

        for cut in self.cutoff_values:           
            mask = self.density_values > cut
            self.voxel_coor = self.xyz_map[mask]             
            self.rho_values = self.density_values[mask]
            self.voxel_indices = self.index_voxels[mask]
            
            # Get the number of clusters(sites) for the new density cutoff
            self.cut_labels = self._clustering_pbc(self.voxel_coor)
            self.n_clusters_cut = len(np.unique(self.cut_labels)) - (1 if -1 in self.cut_labels else 0) 

            if self.n_clusters_cut >= self.n_atoms:                                                 
                self.cluster_spatial_corr = self._sites_cova_pca(self.voxel_coor, self.cut_labels)
              
                max_variance = self.dl/2
                delt_variance = []
                
                for i in range(self.cluster_spatial_corr.shape[1]):
                    cluster_spatial = self.cluster_spatial_corr[:, i].reshape(-1,1)
                    
                    nn1d = NearestNeighbors(n_neighbors = 2).fit(cluster_spatial)
                    distances_nn1d, _ = nn1d.kneighbors(cluster_spatial)
                    distances_nn1d = np.sort(distances_nn1d[:, -1])[-1]
                    delt_variance.append(distances_nn1d)

                if (delt_variance < max_variance).all():
                    # Predict the correct number of type sites
                    all_clusters = np.arange(2, 10)
                    silhouette_averages = [] 
                    labels_shapes_aux =[]

                    if self.scaler == 'MinMaxScaler':
                        data = MinMaxScaler().fit_transform(self.cluster_spatial_corr)  
                    elif self.scaler == 'StandardScaler':
                        data = StandardScaler().fit_transform(self.cluster_spatial_corr)
                    else:
                        sys.exit('Error: Scaler method not recognized') 

                    for n_clusters in all_clusters:                        
                        clustering = KMeans(n_clusters=n_clusters,**kmeans_kwargs)                                             
                        labels = clustering.fit_predict(data) 
                        labels_shapes_aux.append(labels)
                        silhouette_averages.append(silhouette_score(data, labels))                   

                    n_clusters = all_clusters[np.argmax(silhouette_averages)]
                    if (np.max(silhouette_averages) < 0.4):
                        n_clusters = 1
          
                    if n_clusters > 1:
                        self.labels_shapes = labels_shapes_aux[n_clusters - 2]
                    else:
                        self.labels_shapes = np.zeros(len(labels_shapes_aux[0]))

                    upper_size_limit = int(np.ceil(self.maxv/(self.voxel_vol))) # volume of sphere with r=0.90 
                    lower_size_limit = int(np.ceil(self.minv/(self.voxel_vol))) # volume of sphere with r=0.30                    

                    number_voxels =  Counter(self.cut_labels)

                    # Remove the key-value pair with key -1
                    if -1 in number_voxels:
                        del number_voxels[-1]

                    wrong_sites = False
                    # Iterate through the Counter and check values
                    values_arr = []
                    for key, value in number_voxels.items():
                        values_arr.append(value)       
                        if value > upper_size_limit or value < lower_size_limit:                            
                            wrong_sites = True
                            break
                    
                    coordinates = self.voxel_coor[self.cut_labels > -1]
                    labels = self.cut_labels[self.cut_labels > -1]
                                        
                    count_labels_shapes = Counter(self.labels_shapes)
                    # ions% in each of type site could be different?
                    unique_elements = [i for i, count in count_labels_shapes.items() if count <= np.ceil(self.n_atoms * 0.10)]
 
                    if  wrong_sites or len(unique_elements) != 0: 
                        continue
                    else:                  
                        print('Final Density cutoff: {:.4f}'.format(cut)) 
                        self._plot_cluster_shapes(self.cluster_spatial_corr, 'site_shapes_variance.png' ) 

                        smaller_site = 1
                        if smaller_site == 1:
                            self.jump_cut = self.cutoff_values[np.where(self.cutoff_values == cut)[0] + 1]
                        else:
                            self.jump_cut = self.cutoff_values[np.where(self.cutoff_values == cut)[0] + 8]

                        # Write cluster map to be opened with OVITO 
                        coordinates = self.voxel_coor[self.cut_labels > -1]
                        labels = self.cut_labels[self.cut_labels > -1]
                        self._write_cluster_map(coordinates, labels, 'SitesMap_OVITO.dat')              
                        break  
    
        if cut >= value_s and(wrong_sites or len(unique_elements) != 0):
            raise ValueError('\n!!Density cutoff not found!!\n Adjust voxel size or review trajectory\n')

    def _center_sites(self, voxel_coor, labels_cut):     
        max_label = np.max(labels_cut)
        centers = []
        cluster_labels = []
        max_distances = []

        for label in range(max_label + 1):
            cluster_coor = voxel_coor[labels_cut == label]

            for i in range(len(cluster_coor)):
                for j in range(i + 1, len(cluster_coor)):
                    distance = cluster_coor[j] - cluster_coor[i]
                    distance_pos = (distance > 0.5 * self.pbc) * self.pbc
                    distance_neg = (distance < -0.5 * self.pbc) * -self.pbc
                    cluster_coor[j] -= (distance_pos + distance_neg)

            center = np.mean(cluster_coor, axis=0)
            max_distance = np.max(np.linalg.norm(cluster_coor - center, axis = 1))
            
            centers.append(center * self.dl)
            cluster_labels.append(label)
            max_distances.append(max_distance * self.dl[0])        
    
        return np.array(centers), np.array(cluster_labels), np.array(max_distances)

    def _rho_sites(self, m_rho, label_clust):      
        total_rho_cluster = [
            np.sum(m_rho[label_clust == i])
            for i in np.unique(label_clust) if i != -1]

        return np.array(total_rho_cluster)

    def _number_neighbors(self, centers_coor, labels_type):
        
        box_size = self.pbc*self.dl
        centers_coor = np.mod(centers_coor, box_size)
        tree = cKDTree(centers_coor, boxsize=box_size)
        distances, _ = tree.query(centers_coor, k=2)
                
        cutoff_nn = []
        for ty in np.unique(labels_type):
            cutoff =  np.mean(distances[np.where(labels_type == ty)[0], 1]) * 1.20
            cutoff_nn.append(cutoff)        

        nn_neighbors = []
        for id_center in range(len(centers_coor)):
            ty = labels_type[id_center]            
            nn_site = len(tree.query_ball_point(centers_coor[id_center], cutoff_nn[int(ty)])) - 1
            nn_neighbors.append(nn_site)
      
        return cutoff_nn, np.array(nn_neighbors)

    def _write_outputs(self, out_file, info):        
        header = '#Site-label   Center(X)   Center(Y)   Center(Z)  Site-Type   N-Neighbors   Occupancy  Site-Amplitud'

        with open(out_file, 'w') as file:
            file.write(header + "\n")
            i = 1
            for row in info:                
                formatted_row = '{:^12}{:^12.6f}{:^12.6f}{:^12.6f}{:^12}{:^12}{:^13.6f}{:^13.6f}'.format(
                    i, float(row[0]), float(row[1]), float(row[2]),
                    int(row[3]), int(row[4]), float(row[5]), float(row[6])                    
                )
                file.write(formatted_row + "\n")
                i += 1

    def _read_voxel_index(self, voxels_visited_file):        
        voxel_indices = []

        with open(voxels_visited_file, 'r') as file:
            next(file)  # Skip the first line (header)
            for line in file:               
                values = [int(val) for val in line.split()]
                voxel_indices.append(values)

        return np.array(voxel_indices)
    
    def _atom_jumps(self, voxel_indices_time, cut_labels, voxel_indices,\
                     n_frames):
        
        def interstice_remover(arr):
            # Number of rows and columns
            rows = len(arr)
            cols = len(arr[0])

            # Iterate over each row
            for row in range(rows):
                for col in range(cols):
                    if arr[row][col] == -1:
                        # Try to replace with a suitable previous value
                        replaced = False

                        # Check previous values in the row
                        for prev_col in range(col-1, -1, -1):
                            if arr[row][prev_col] >= 0 and all(arr[r][col] != arr[row][prev_col] for r in range(rows)):
                                arr[row][col] = arr[row][prev_col]
                                replaced = True
                                break

                        # If no replacement found, try next values
                        if not replaced:
                            for next_col in range(col+1, cols):
                                if arr[row][next_col] >= 0 and all(arr[r][col] != arr[row][next_col] for r in range(rows)):
                                    arr[row][col] = arr[row][next_col]
                                    break
            return arr       

        def find_connected_jumps(graph, simultaneous_jumps):
            visited = set()
            components = []

            for node in graph:
                if node not in visited:
                    stack = [node]
                    component = set()
                    component_pairs = []

                    while stack:
                        current = stack.pop()
                        if current not in visited:
                            visited.add(current)
                            component.add(current)
                            stack.extend(graph[current] - visited)
                            
                    for pair in simultaneous_jumps:
                        if pair[0] in component or pair[1] in component:
                            component_pairs.append(pair)
                    components.append(component_pairs)

            return components

        def find_nearest_non_negative_left(arr, index):
            left = index
            while left >= 0:
                if arr[left] != -1:
                    return arr[left]
                left -= 1
            return None

        def find_nearest_non_negative_right(arr, index):
            right = index
            while right < len(arr):
                if arr[right] != -1:
                    return arr[right]
                right += 1
            return None
              
        mask = self.density_values > self.jump_cut
        voxel_coor = self.xyz_map[mask]        
        voxel_indices = self.index_voxels[mask]
        cut_labels = self._clustering_pbc(voxel_coor)
        n_clusters_cut = len(np.unique(cut_labels)) - (1 if -1 in cut_labels else 0) 

        if n_clusters_cut != self.n_clusters_cut:
            sys.exit(f"Error: Site count mismatch. Jump analysis: {n_clusters_cut}, Sites analysis: {self.n_clusters_cut}.")

        atoms_jumps_indices = []; atoms_jumps_reverse_indices = []
        njumps_time = []; njumps_time_noreverse = []
        sites_visited_atom = []; atoms_paths = []; omega = [] 
        index_label_map = {index: label for index, label in zip(voxel_indices, cut_labels)}
                
        nframes_allow = (n_frames // self.deltf) * self.deltf     
        
        for indices_atom in voxel_indices_time:
            mapped_labels_time = [index_label_map.get(index, -1) for index in indices_atom]
            atom_path = np.array(mapped_labels_time)
            atom_path = atom_path[:nframes_allow]

            if self.deltf > 1:
                atom_path_reshape = atom_path.reshape(-1, self.deltf)
                atom_path = atom_path_reshape[:, 0] 

            atoms_paths.append(atom_path)

        atoms_paths_filtered = interstice_remover(atoms_paths)   

        # Calculate how many frames have atoms in interstice (-1 values)
        negative_counts = [np.sum(row == -1) for row in atoms_paths_filtered if -1 in row]
        if negative_counts:
            average_negatives_per_row = np.mean(negative_counts)
        else:
            average_negatives_per_row = 0 

        if average_negatives_per_row > 10:         
            sys.exit(f'Suggestion: Increase -deltf (current -1 average count per line: {average_negatives_per_row})')   

        for atom_path in atoms_paths_filtered:            
            # Calculating simultaneous jumps
            atom_paths_not_noise = atom_path[atom_path >= 0]
            atom_paths_original_indices = np.where(atom_path >= 0)[0]  
            site_change = np.where(atom_paths_not_noise[:-1] != atom_paths_not_noise[1:])[0]
            
            sites_changes_i = atom_paths_original_indices[site_change]
            sites_changes_f = atom_paths_original_indices[site_change + 1]

            # simultaneous jumps calculation ends
            jump_site_labels_i = atom_path[sites_changes_i]
            jump_site_labels_f = atom_path[sites_changes_f]
            
            sites_changes_back = np.where(jump_site_labels_i[:-1] == jump_site_labels_f[1:])[0]   
            
            jumps_arr = np.zeros_like(atom_path)
            for start, end in zip(sites_changes_i , sites_changes_f):
                if end == start + 1:
                    jumps_arr[start] = 1
                jumps_arr[start + 1:end] = 1

            omega.append(atom_path)

            njumps_time.append(jumps_arr)
            njumps_time_noreverse.append(jumps_arr)

            atoms_jumps_indices.append(len(np.unique(site_change)))            
            atoms_jumps_reverse_indices.append(len(np.unique(sites_changes_back)))

            sites_visited_atom.append(len(np.unique(atom_paths_not_noise)))

        #Number of simultaneous jumps according "prl 116, 135901 (2016)".
        ns_jumps = np.sum(np.array(njumps_time), axis=0)
  
        # Calculate strings of ions according to "prl 116, 135901 (2016)".
        sites_atoms_time = np.array(omega) 
    
        self.total_time_updated = ((nframes_allow-1) // self.deltf) * (self.time_frames * self.deltf) 
    
        #Matrix of jumps. Fix of last frame jump.
        njumps_time_noreverse = np.array(njumps_time_noreverse)
  
        traj_chunks = sites_atoms_time
      
        traj_jumps = njumps_time_noreverse 
   
        sites_atoms_time = traj_chunks
        jumps_atoms_time = traj_jumps
        jumps_atoms_time[:,-1] = 0

        forwardjumps = np.sum(np.array(jumps_atoms_time), axis=0)
        string_indices_time = np.where(forwardjumps > 0)[0]

        # Find the indices of the rows for which the value is 1 in new_arr at arr_jumps_index
        rows_with_jumps = {}
        for col in string_indices_time:
            rows_with_jumps[col] = [i for i, row in enumerate(jumps_atoms_time) if row[col] == 1]
        
        string_lengths = []; long_strings = []
        for col, rows in rows_with_jumps.items():
            #print(f"\nTime step {col}:")
            simultaneous_jumps_time = []
            for row in rows:                 
                element_col = sites_atoms_time[row, col]
                element_col_plus_1 = sites_atoms_time[row, col + 1]
                if element_col == -1 or element_col_plus_1 == -1:
                    left_value = find_nearest_non_negative_left(sites_atoms_time[row], col)
                    right_value = find_nearest_non_negative_right(sites_atoms_time[row], col + 1)
                    jump_sites = [left_value, right_value]                    
                else:
                    jump_sites = [element_col, element_col_plus_1]
                
                simultaneous_jumps_time.append(jump_sites)
                            
            graph = defaultdict(set)
            #Build the graph
            for a, b in simultaneous_jumps_time:
                graph[a].add(b)
                graph[b].add(a)                   

            connected_jumps = find_connected_jumps(graph, simultaneous_jumps_time)
            
            cluster_sizes = [len(cluster) for cluster in connected_jumps]
            
            cluster_size_counts = Counter(cluster_sizes)
        
            string_lengths.append(cluster_size_counts)

        strings_per_length = sum((Counter(counter) for counter in string_lengths), Counter())

        prob_string = {key: ((value)/ (self.total_time_updated))for key, value in strings_per_length.items()}
    
        return  np.array(atoms_jumps_indices),\
                np.array(atoms_jumps_reverse_indices), ns_jumps,\
                prob_string, sites_visited_atom
 
    def _write_string_freq(self, prob_string, output_filename):
        
        sorted_items = sorted(prob_string.items())     
        
        with open(output_filename, 'w') as file:
            file.write("# String Length (#n)  f_n(ps^-1)\n")
            for key, value in sorted_items:
                file.write(f"{key:10.4f}  {value:10.4f}\n")
    
    def _write_jumps(self, Njumps, RNjumps, SNjumps, output_jumpsinfo_filename, output_stringjumps_filename): 
        self_corr = np.divide(RNjumps, Njumps, where=Njumps!=0, out=np.full_like(RNjumps, np.nan, dtype=float))
        self_corr = np.nan_to_num(self_corr, nan=-1)  

        with open(output_jumpsinfo_filename, 'w') as file:
            file.write("# Atom      NJumps     RNJumps    JumpCoefficient\n")
            i = 1
            for row in range(len(Njumps)):
                formatted_row = '{:^12}{:^12}{:^12}{:^12.4f}'.format(
                    i, int(Njumps[row]), int(RNjumps[row]), float(self_corr[row])
                )
                file.write(formatted_row + "\n")
                i += 1
        
        with open(output_stringjumps_filename, 'w') as file:
            file.write("# Frame     NSjumps \n")
            i = 1
            for row in range(len(SNjumps)):
                formatted_row = '{:^12}{:^12}'.format(
                    i, SNjumps[row]
                )
                file.write(formatted_row + "\n")
                i += 1
               
    def _plot_cluster_data(self, amplitud_values, output_filename):
        
        amplitud_indices = range(len(amplitud_values))
        
        fig, ax = plt.subplots(figsize=(8.0, 8.0))        
        ax.scatter(amplitud_indices, amplitud_values, c = self.labels_shapes, cmap='rainbow', edgecolor='k', s = 40)
        ax.set_xlabel('Site Index', fontproperties=self.font_properties, fontsize=30)
        ax.set_ylabel('Site amplitude (A)', fontproperties=self.font_properties, fontsize=30)
        ax.tick_params(axis='both', labelsize=20, labelrotation=0)

        fig.tight_layout()        
        plt.savefig(output_filename)
        plt.close()
    
    def _plot_cluster_shapes(self, variances, output_filename):
        title = 'Site shapes variance'
        fig = plt.figure(figsize=(10.0, 8.0))
        ax = fig.add_subplot(111, projection='3d')   

        ax.scatter(variances[:, 0],
                   variances[:, 1],
                   variances[:, 2],
                   c = self.labels_shapes,
                   cmap='rainbow', edgecolor='k', s=40)
            
        ax.set_title(title,fontproperties=self.font_properties, fontsize=20, fontweight='bold')
        ax.set_xlabel('x', fontproperties=self.font_properties, fontsize=14)
        ax.set_ylabel('y', fontproperties=self.font_properties, fontsize=14)
        ax.set_zlabel('z', fontproperties=self.font_properties, fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=12)
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

        fig.tight_layout()        
        plt.savefig(output_filename, bbox_inches='tight')
        plt.close()

    def _plot_sites_visited(self, sites_visited):
        #Ploting number of sites visited per atom
        num_bins = int(np.sqrt(len(sites_visited)))
        fig, ax = plt.subplots(figsize=(8.0, 8.0))   

        sites_visited_corr = np.array(sites_visited) 

        ax.hist(sites_visited_corr,
                bins = num_bins,
                color = 'b', alpha=0.7, edgecolor='k')
        ax.set_xlabel('Number of sites visited', fontproperties = self.font_properties, fontsize=30) 
        ax.set_ylabel('Number of Atoms',fontproperties = self.font_properties, fontsize=30)

        # Ensure y-axis has integer values only
        ax.xaxis.set_major_locator(MaxNLocator(integer=True)) 
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))        

        ax.tick_params(axis='both', labelsize=20, labelrotation=0)
        fig.tight_layout()        
        plt.savefig('sites_visited.png')
        plt.close()
        
    def _process_extra_data(self, density_maps, sites_outputs, voxels_visited, jumpsinfos,
                            simultaneousjumps, stringfreq):     
        for density_map, sites_output, voxel_path, jumpsinfo, simultaneousjump, stringfreq in\
              zip(density_maps, sites_outputs, voxels_visited, jumpsinfos, simultaneousjumps, stringfreq):
            print('\nReading file:', density_map)
            density_values, _ , index_voxels = self._read_density_map(density_map)
            
            # Copying labels and indices from the complete density map                    
            mask_indices_ref = np.isin(self.voxel_indices, index_voxels)
            mask_indices_tar = np.isin(index_voxels, self.voxel_indices)

            voxel_indices = self.voxel_indices[mask_indices_ref]
            cut_labels = self.cut_labels[mask_indices_ref]  
            cluster_labels = np.unique(cut_labels)

            label_shape = self.labels_shapes[np.isin(self.site_labels, cluster_labels)]
            rho_site = self._rho_sites(density_values[mask_indices_tar], cut_labels)
            center_sites = self.center_sites[np.isin(self.site_labels, cut_labels)]
            amplitud_sites = self.amplitud_sites[np.isin(self.site_labels, cut_labels)]

            _ , nn_neighbors = self._number_neighbors(center_sites, label_shape)
           

            info = np.column_stack([
                center_sites[:,0], center_sites[:,1], center_sites[:,2],              
                label_shape, 
                nn_neighbors,              
                rho_site * self.voxel_vol, 
                amplitud_sites
            ])

            print('Writing file:', sites_output)
            self._write_outputs(sites_output, info)

            # Analyzing and printing information per site type            
            print(f'\nInformation for site type')
            for site_type in np.unique(info[:, 3]):
                site_specific_info = info[info[:, 3] == site_type]
                print(f'Label site: {int(site_type)}, Number of Sites: {len(site_specific_info)}')

            # Jump imformation for cluster
            print(f'\nReading file: {voxel_path}')
            voxel_indices_time = self._read_voxel_index(voxel_path)
            print(f'Total time: {self.time_frames * (voxel_indices_time.shape[1]-1):.4f}')
            print(f'Time step VoxelIndices: {self.time_frames}')

            total_time = (voxel_indices_time.shape[1]-1) * (self.time_frames)
            total_jupms, reverse_jumps, simultaneous_jumps, average_probabilities, _=\
                  self._atom_jumps(voxel_indices_time, cut_labels, voxel_indices,\
                                   voxel_indices_time.shape[1], total_time)            
            
            print(f'\nWriting file: {jumpsinfo}')

            print(f'\nWriting file: {simultaneousjump}')
            print(f'Frames interval: {self.deltf:d}')
            print(f'Time interval between frame intervals: {self.time_frames * self.deltf:.4f}')
            print(f'Updated total time: {self.total_time_updated:.4f}')
            
            self._write_jumps(total_jupms, reverse_jumps, simultaneous_jumps,\
                               jumpsinfo, simultaneousjump)


            print(f'\nWriting file: {stringfreq}')                  
            self._write_string_freq(average_probabilities, stringfreq)

    def _site_analysis(self):
        try:
            self.voxels_visited_file = 'VoxelIndices.dat'            
            self.voxel_indices_time = self._read_voxel_index(self.voxels_visited_file)           
            self.n_atoms = self.voxel_indices_time.shape[0]
            self.n_frames = self.voxel_indices_time.shape[1]
          
            self._find_density_cutoff()
            # Calculate the center of each cluster and the amplitud         
            self.center_sites, self.site_labels, self.amplitud_sites = self._center_sites(self.voxel_coor, self.cut_labels)
            
            print('\nAmplitude types: {:d}'.format(len(np.unique(self.labels_shapes))))

            am_values = []

            for am in range(len(np.unique(self.labels_shapes))):
                mask = self.labels_shapes == am
                am_avg = np.mean(self.amplitud_sites[mask])               
                am_value = am_avg 
                am_values.append(am_value)
                print(f'Site amplitudes type {am}: {am_value:.4f}')

            self.cutoff_nneighbors, self.nn_neighbors = self._number_neighbors(self.center_sites, self.labels_shapes)
            
            for cut_nn in range(len(np.unique(self.labels_shapes))):
                cutoff_nn_value = self.cutoff_nneighbors[cut_nn]             
                print(f'Distance used as cutoff in NN type {cut_nn}: {cutoff_nn_value:.4f}')
            
            self.rho_site = self._rho_sites(self.rho_values, self.cut_labels)
            
            # Preparing data
            info = np.column_stack([
                self.center_sites[:,0], self.center_sites[:,1], self.center_sites[:,2],              
                self.labels_shapes, 
                self.nn_neighbors,              
                self.rho_site * self.voxel_vol, 
                self.amplitud_sites,             
            ])
                
            # Writing information file
            print('\n>> Sites map of the structure <<')
            print(f'Writing file: {self.out_file}')
            self._write_outputs(self.out_file, info)

            # Analyzing and printing information per site type          
            print('\nInformation for site type')
            for site_type in np.unique(info[:, 3]):
                site_specific_info = info[info[:, 3] == site_type]              
                print(f'Label site: {int(site_type)}, Number of Sites: {len(site_specific_info)}')
            
            self._plot_cluster_data(self.amplitud_sites, 'site_amplitudes.png')

            self.total_jupms, self.reverse_jumps, self.simultaneous_jumps,\
                  self.average_probabilities, self.sites_visited_atom =\
                  self._atom_jumps(self.voxel_indices_time, self.cut_labels,\
                                    self.voxel_indices, self.n_frames)

            self._plot_sites_visited(self.sites_visited_atom)

            print(f'\nReading file: {self.voxels_visited_file}')
            print('>> Trajectory data <<')
            print(f'Total time: {self.total_time:.4f}')            
            print(f'Time step VoxelIndices: {self.time_frames}')
                       
            ############   
            self.jumps_file = 'jumps_info.dat'            
            self.simultaneous_jumps_file = 'simultaneous_jumps.dat'
            self.string_freq_file = 'string_frequency.dat'

            print(f'\nWriting file: {self.jumps_file}')
            print(f'Writing file: {self.simultaneous_jumps_file}')
            print(f'Writing file: {self.string_freq_file}') 
            

            print(f'\nFrames interval: {self.deltf:d}')
            print(f'Time interval between frame intervals: {self.time_frames * self.deltf:.4f}')            
            print(f'Updated total time: {self.total_time_updated:.4f}')

            self._write_jumps(self.total_jupms, self.reverse_jumps, self.simultaneous_jumps,\
                               self.jumps_file, self.simultaneous_jumps_file)
                   
            self._write_string_freq(self.average_probabilities, self.string_freq_file)
          
            if self.clustering == 1 and self.doped == 0:
                print('\n>> Sites map of <di> clusters <<')

                cluster_labels = ['H-', 'M-', 'L-', 'HM-'] if self.verbosity == 1 else ['HM-', 'L-']
                density_maps_outputs = [f'{suffix}{self.density_map_targ}' for suffix in cluster_labels]           
                sites_outputs = [f'{suffix}{self.out_file}' for suffix in cluster_labels]
                voxels_visited_outputs = [f'{suffix}{self.voxels_visited_file}' for suffix in cluster_labels]
                jumpsinfo_outputs = [f'{suffix}{self.jumps_file}' for suffix in cluster_labels]
                simultaneousjumps_outputs = [f'{suffix}{self.simultaneous_jumps_file}' for suffix in cluster_labels]
                stringfreq_outputs = [f'{suffix}{self.string_freq_file}' for suffix in cluster_labels]

                self._process_extra_data(density_maps_outputs, sites_outputs, voxels_visited_outputs,
                                          jumpsinfo_outputs, simultaneousjumps_outputs, stringfreq_outputs)

            elif self.clustering == 0 and self.doped == 1:
                print('\n>> Sites map of the doped structure <<')
                
                cluster_labels = ['Dope-']
                density_maps_outputs = [f'{suffix}{self.density_map_targ}' for suffix in cluster_labels]           
                sites_outputs = [f'{suffix}{self.out_file}' for suffix in cluster_labels]
                voxels_visited_outputs = [f'{suffix}{self.voxels_visited_file}' for suffix in cluster_labels]
                jumpsinfo_outputs = [f'{suffix}{self.jumps_file}' for suffix in cluster_labels]
                simultaneousjumps_outputs = [f'{suffix}{self.simultaneous_jumps_file}' for suffix in cluster_labels]
                stringfreq_outputs = [f'{suffix}{self.string_freq_file}' for suffix in cluster_labels]

                self._process_extra_data(density_maps_outputs, sites_outputs, voxels_visited_outputs,
                                          jumpsinfo_outputs, simultaneousjumps_outputs, stringfreq_outputs)
                
            elif self.clustering == 1 and self.doped == 1:
                print('\n>> Sites map of the doped structure and its <di> clusters <<')
                
                cluster_labels = ['Dope-','H-Dope-', 'M-Dope-', 'L-Dope-', 'HM-Dope-'] if self.verbosity == 1 else ['Dope-','HM-Dope-', 'L-Dope-']
                density_maps_outputs = [f'{suffix}{self.density_map_targ}' for suffix in cluster_labels]           
                sites_outputs = [f'{suffix}{self.out_file}' for suffix in cluster_labels]
                voxels_visited_outputs = [f'{suffix}{self.voxels_visited_file}' for suffix in cluster_labels]
                jumpsinfo_outputs = [f'{suffix}{self.jumps_file}' for suffix in cluster_labels]
                simultaneousjumps_outputs = [f'{suffix}{self.simultaneous_jumps_file}' for suffix in cluster_labels]
                stringfreq_outputs = [f'{suffix}{self.string_freq_file}' for suffix in cluster_labels]

                self._process_extra_data(density_maps_outputs, sites_outputs, voxels_visited_outputs,
                                          jumpsinfo_outputs, simultaneousjumps_outputs, stringfreq_outputs)

        except ValueError as e:
            print(e)
