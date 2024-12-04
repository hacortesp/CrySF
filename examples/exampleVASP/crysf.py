import argparse
from libraries.CrySF import CrySF 
from timeit import default_timer as timer

def parse_args():
    parser = argparse.ArgumentParser(description='Parse input arguments for the data processing script.')    
    parser.add_argument('-nts'    ,  type = float, help='Time interval between frames in VoxelIndices.dat (-ts * -tts)')
    parser.add_argument('-minv'   ,  type = float, default = 0.28, help='Minimum volume condition for the site.')    
    parser.add_argument('-maxv'   ,  type = float, default = 3.05, help='Maximum volume condition for the site.')      
    parser.add_argument('-verb'   ,  type = int,   choices = [0, 1], default=0, help='Set verbosity level, 0 = L and HM(together) or 1 = L, M, H and HM(together).')
    parser.add_argument('-deltf'  ,  type = int,   default = 1, help='Eevery how many frames the trajectory is used to calculate the jump coeficient, simultaneous jumps.')
    parser.add_argument('-scaler' ,  type = str,   default = 'MinMaxScaler', help='Scaler method used for site type identification. Options: MinMaxScaler, StandardScaler.')
    parser.add_argument('-clus'   ,  type = int,   choices = [0, 1], default=0, help='Set whether clustering <di> of the trajectory is applied (0 = NO or 1 = YES).')
    parser.add_argument('-dop'    ,  type = int,   choices = [0, 1], default=0, help='Sets whether to use a doped density map, 0 = NO or 1 = YES.')
    return parser.parse_args()

args = parse_args()
start = timer()

CrySF(
   args.nts, args.minv, args.maxv, args.verb, args.deltf, args.scaler, args.clus, args.dop
)

end = timer()
elapsed_time = (end - start) / 60
print(f"\ntime elapsed = {elapsed_time:.4f} minutes")

