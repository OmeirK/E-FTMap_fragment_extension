'''
Written by: Omeir Khan

Split 2d screening hits into different files based on their heavy atom
count (HAC). Calculate number of rotatable bonds in each molecule, and 
store these values in the output file
'''

import os
import argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

parser = argparse.ArgumentParser()

parser.add_argument('--screening_hits', '-s', help='The output 2d_stage_1.txt file with screening hits, sorted by SMARTS coverage and tab separated.')
parser.add_argument('--outdir', '-od', help='The name of a directory where sorted 2d hits will be stored. If the folder does not exist it will be created. (default = sorted_2d_hits/)', default='sorted_2d_hits/')

args = parser.parse_args()

def read_screening_hits(hit_list):
    output_data = {}
    with open(hit_list) as f:
        for line in f:
            l_data = line.strip().split('\t')
            smi = l_data[0]
            name = l_data[1]
            smarts = l_data[2]

            mol = Chem.MolFromSmiles(smi)
            
            if mol is None:
                continue

            hac = mol.GetNumHeavyAtoms()
            rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

            if hac not in output_data:
                outfile = f'{args.outdir}/2d_hits_{hac}_HAC.tsv'
                if os.path.exists(outfile):
                    os.remove(outfile)

                output_data[hac] = outfile
            
            with open(output_data[hac], 'a') as fo:
                fo.write(f'{smi}\t{name}\t{smarts}\t{hac}\t{rot_bonds}\n')

def main():
    os.makedirs(f'{args.outdir}', exist_ok=True)
    read_screening_hits(args.screening_hits)

if __name__=='__main__':
    main()
