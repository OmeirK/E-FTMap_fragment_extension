'''
Collect the subset of 2d hits that match the specified criteria
'''

import os
import glob
import argparse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

parser = argparse.ArgumentParser()

parser.add_argument('--query_mol', '-m', help='The .mol file with the aligned query mol ')
parser.add_argument('--sorted_hits_dir', '-sd', help='The directory containing hits sorted by HAC')
parser.add_argument('--outfile', '-o', help='A .tsv file with the filtered hits (default = filtered_hits.tsv)', default='filtered_hits.tsv')
parser.add_argument('--n_new_rot_bonds', '-nb', help='The maximum number of additional rotatable bonds that can be included in filtered compounds. The baseline number of rotatble bonds is determined by the --query_mol that is provided. (default = 2)', default=2, type=int)
parser.add_argument('--max_hac', '-maxhac', help='The maximum number of heavy atoms that can be included in filtered compounds (default = 38)', default=38, type=int)
parser.add_argument('--min_hac', '-minhac', help='The minimum number of heavy atoms that can be included in filtered compounds (default = 0)', default=0, type=int)

args = parser.parse_args()

def get_subset(hit_f, n_new_rot_bonds, n_frag_rot_bonds, min_hac, max_hac):
    hits = []
    bond_lim = n_frag_rot_bonds + n_new_rot_bonds

    with open(hit_f) as f:
        for line in f:
            l_data = line.strip().split('\t')
            n_rot = int(l_data[-1])
            hac = int(l_data[-2])

            if n_rot > bond_lim:
                continue
            else:
                if (hac >= min_hac) and (hac <= max_hac):
                    hits.append(line)
    return hits

def main():
    frag_mol = Chem.MolFromMolFile(args.query_mol)
    frag_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(frag_mol)
    frag_hac = frag_mol.GetNumHeavyAtoms()
    print(frag_rot_bonds, frag_hac)

    sorted_hits = glob.glob(f'{args.sorted_hits_dir}/*.tsv')
    sorted_hits.sort()

    print(sorted_hits)
    
    filtered_hits = []
    for hit_f in sorted_hits: 
        hit_f_name = os.path.basename(hit_f)
        hit_f_hac = int(hit_f_name.split('_')[2])

        print(hit_f_hac, hit_f_name)

        hits = get_subset(hit_f, args.n_new_rot_bonds, frag_rot_bonds, args.min_hac, args.max_hac)

        filtered_hits += hits

        print(hit_f_name, len(hits))

    with open(args.outfile, 'w') as fo:
        fo.write(''.join(filtered_hits))



if __name__=='__main__':
    main()
