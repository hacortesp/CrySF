# main.py
import argparse
from libraries.CrySF import CrySF
from timeit import default_timer as timer

def parse_args():
    parser = argparse.ArgumentParser(description='Parse input arguments for the data processing script.')    
    parser.add_argument('-nts'    , type=float, help='Time interval between frames in VoxelIndices.dat (-ts * -tts)')
    parser.add_argument('-minv'   , type=float, default=0.28, help='Minimum volume condition for the site.')
    parser.add_argument('-maxv'   , type=float, default=3.05, help='Maximum volume condition for the site.')
    parser.add_argument('-verb'   , type=int, choices=[0, 1], default=0, help='Verbosity: 0 = L & HM or 1 = L, M, H & HM.')
    parser.add_argument('-deltf'  , type=int, default=1, help='Every how many frames to compute jump coefficient.')
    parser.add_argument('-struct' , type=str, choices=['pure', 'defect'], default='pure', help='Structure type: pure or defect.')
    parser.add_argument('-scaler' , type=str, default='MinMaxScaler', help='Scaler for site type identification.')
    parser.add_argument('-clus'   , type=int, choices=[0, 1], default=0, help='Clustering <di> in trajectory: 0 = NO, 1 = YES.')
    parser.add_argument('-dop'    , type=int, choices=[0, 1], default=0, help='Use doped density map: 0 = NO, 1 = YES.')
    return parser.parse_args()

args = parse_args()
start = timer()

CrySF(
    args.nts,
    args.minv,
    args.maxv,
    args.verb,
    args.deltf,
    args.struct,
    args.scaler,
    args.clus,
    args.dop
)

end = timer()
elapsed_time = (end - start) / 60
print(f"\ntime elapsed = {elapsed_time:.4f} minutes")
