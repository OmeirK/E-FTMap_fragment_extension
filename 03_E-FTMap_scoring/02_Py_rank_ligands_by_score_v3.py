'''
Written By: Omeir Khan

Once ligands have been scored, rank them by score. 
Save the best scoring pose from each ligand

v2: [X] Add an argument to save the top N ligands
v3: [X] Consider scores from multiple tsv files. To do
       this it must map score files to sdf files
    
'''

import os
import argparse
import pandas as pd
from rdkit import Chem

parser = argparse.ArgumentParser()

#parser.add_argument('--eftmap_score_tsv', '-s', help='A .tsv file containing E-FTMap scores for compounds')
parser.add_argument('--eftmap_score_tsv_dir', '-sd', help='A directory containing .tsv files with E-FTMap scores for compounds')
parser.add_argument('--out_prefix', '-op', help='A prefix to add before the names of output files (default = ranked)', default='ranked')
parser.add_argument('--n_ligands', '-n', help='Specify the number of top ranked ligands to save to an SDF file (default = 100)', default=100, type=int)

args = parser.parse_args()

# Save data for the top scoring pose for each ligand.
# This is based on the ligand names! If they are not
# provided in the original file this wont work!
def read_tsv(tsv_file, data_dict = {}):
    read_sdfs = {}

    with open(tsv_file) as f:
        tsv_lines = f.readlines()

    header = tsv_lines[0]

    df = pd.read_csv(tsv_file, delimiter='\t')
    for i, ligname in enumerate(df['Lig_Name']):
        #print(i, ligname)
        if ligname not in data_dict:
            data_dict[ligname] = {'score': -1, 'mol_idx': -1, 'line': ''}
        
        origin_sdf = df['Origin_File'][i]
        mol_idx = int(df['Mol_Idx'][i])
        score = float(df['Total'][i])
        
        if origin_sdf not in read_sdfs:
            mols = Chem.SDMolSupplier(origin_sdf)
            read_sdfs[origin_sdf] = mols
        
        mol = read_sdfs[origin_sdf][mol_idx]

        if score > data_dict[ligname]['score']:
            data_dict[ligname]['score'] = score
            data_dict[ligname]['mol_idx'] = mol_idx
            data_dict[ligname]['mol'] = mol
            data_dict[ligname]['line'] = tsv_lines[i+1]
    
    return data_dict, header

def main():
    # Save sorted data as an sdf file
    score_data = {}
    tsv_files = os.listdir(args.eftmap_score_tsv_dir)
    for tsv in tsv_files:
        tsv_path = f'{args.eftmap_score_tsv_dir}/{tsv}'
        print(f'Reading {tsv_path}...')
        score_data, header = read_tsv(tsv_path)


    # Sort the data by score, from greatest to least
    score_data = dict(sorted(score_data.items(), key=lambda k_v: k_v[1]['score'], reverse=True))

    fo = open(f'{args.out_prefix}_score_list.tsv', 'w')
    fo.write(header)

    n = 0 # Log number of ligands saved to an sdf file
    with Chem.SDWriter(f'{args.out_prefix}_top{args.n_ligands}_mols.sdf') as w:
        for ligname in score_data:
            fo.write(score_data[ligname]['line'])
            mol_idx = score_data[ligname]['mol_idx']
            mol = score_data[ligname]['mol']

            if n < args.n_ligands:
                print(ligname, score_data[ligname]['score'])
                w.write(score_data[ligname]['mol'])
                n += 1

    fo.close()


if __name__=='__main__':
    main()
