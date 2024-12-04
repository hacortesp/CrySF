import argparse
from libraries.DensityMap import DensityMap
from timeit import default_timer as timer

def parse_args():
    parser = argparse.ArgumentParser(description='Input data')
    parser.add_argument('-to'  , type=str,   default='topo.gro', help='Name of topology file (.gro or .data formats)')
    parser.add_argument('-tr'  , type=str,   default='trajectory.xtc', help='Name of trajectory file')
    parser.add_argument('-f'   , type=str,   default='XDATCAR', help='Format options: XTC, TRR, LAMMPS, XDATCAR')
    parser.add_argument('-ts'  , type=float, default=0.15, help='Time interval between frames in the trajectory')
    parser.add_argument('-v'   , type=float, default=0.2, help='Voxel size')
    parser.add_argument('-a'   , type=str,   default='Li', help='Atom selection: integer (index) or string (name) based on trajectory format')
    parser.add_argument('-tts' , type=int,   default= 1, help='Frames interval to print the VoxelIndices.dat file')
    parser.add_argument('-verb', type=int,   choices=[0, 1], default= 0, help='Set verbosity level (0 = L and HM or 1 = L, M, H and HM)')
    parser.add_argument('-clus', type=int,   choices=[0, 1], default= 0, help='Set whether clustering <di> of the trajectory is applied (0 = NO or 1 = YES)')
    
    
    return parser.parse_args()


args = parse_args()
start = timer()

DensityMap(args.to, args.tr, args.f, args.ts, args.v, args.a, args.tts, args.verb, args.clus)

end = timer()
elapsed_time = (end - start) / 60
print(f"time elapsed = {elapsed_time:.4f} minutes")

