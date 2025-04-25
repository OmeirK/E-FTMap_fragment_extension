import os
import json
import shutil
import argparse
import subprocess
import multiprocessing as mp

parser = argparse.ArgumentParser()

parser.add_argument('--job_file', '-f', help='The job file to be used for 3D placement')
parser.add_argument('--recpdb', '-r', help='A .pdb files of the aligned receptor to dock to')
parser.add_argument('--ligmol', '-l', help='A .mol file with the coordinates of the template ligand to use for docking')
parser.add_argument('--outdir', '-od', help='Name of an directory where the output will be stored')
parser.add_argument('--n_cpus', '-cpu', help='Number of cpus to use. Increase to parallelize the calculation (default=28)', default=28, type=int)
parser.add_argument('--n_confs', '-nc', help='(Optional) The number of conformers to generate for template-based docking (default=100)', default=100, type=int)

args = parser.parse_args()

CMD_PATH='/projectnb/docking/omeir/E-FTMAP_molecular_stitching/git_repo/E-FTMap_fragment_extension/02_3d_filtering/'

def mp_func(cmd, ligid):
    cmd = cmd.split('\t')
    subprocess.run(cmd)

def main():
    # Create the output directory
    os.makedirs(args.outdir, exist_ok=True)

    
    # Read the job file
    command_list = []
    with open(args.job_file) as f:
        for line in f:
            l_data = line.strip().split('\t')
            lig_smi = l_data[0]
            lig_name = l_data[1]
            smarts = l_data[2]

            cmd = f'python3\t{CMD_PATH}/ligand_3d_placement_v6.py\t-r\t{args.recpdb}\t-tl\t{args.ligmol}\t-as\t{smarts}\t-ls\t{lig_smi}\t-ln\t{lig_name}\t-od\t{args.outdir}\t-nrc\t-nc\t{args.n_confs}'
            print(cmd)
            #command_list.append(cmd.split('\t'))
            command_list.append((cmd, lig_name))

    
    print(len(command_list))
    with mp.Pool(args.n_cpus) as pool:
        pool.starmap(mp_func, command_list, chunksize=1)

if __name__=='__main__':
    main()

