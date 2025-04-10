'''
Written By: Omeir Khan

This script will take all sdf files in a folder, and combine them into one large file.
This is useful for processing the output of the 3d generation stage after the 2d screen, 
as it can be used prior to scoring poses. Note that if there are more than 50,000
poses to score, the data is split across multiple .sdf files to avoid memory issues
'''

import os
import tqdm
import argparse
from rdkit import Chem

parser = argparse.ArgumentParser()

parser.add_argument('--sdf_dir', '-d', help='Directory with sdf files outputted from the 3D placement calculation')
parser.add_argument('--outdir', '-od', help='Name of the output directory (default = combined_poses/)', default='combined_poses/')

args = parser.parse_args()

def main():
    os.makedirs(args.outdir, exist_ok=True)
    #sdf_list = glob.glob(f'{args.sdf_dir}/*.sdf')

    all_mols = []
    
    sdf_list = os.listdir(f'{args.sdf_dir}')
    out_idx = 0
    
    all_mols = []
    for sdf in tqdm.tqdm(sdf_list):
        mols = Chem.SDMolSupplier(f'{args.sdf_dir}/{sdf}')
        all_mols += mols

        if len(all_mols) > 50000:
            subset_outfile = f'{args.outdir}/combined_poses_{out_idx}.sdf'
            print(f'Save {out_idx} subset to output {subset_outfile}')

            with Chem.SDWriter(subset_outfile) as w:
                for m in all_mols:
                    w.write(m)
            
            all_mols = []
            out_idx  += 1

    subset_outfile = f'{args.outdir}/combined_poses_{out_idx}.sdf'
    with Chem.SDWriter(subset_outfile) as w:
        for m in all_mols:
            w.write(m)
    
    all_mols = []
    out_idx  += 1

if __name__=='__main__':
    main()
